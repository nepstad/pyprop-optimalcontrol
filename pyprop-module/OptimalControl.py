import warnings
from scipy import sparse

class OptimalControl:
	"""
	Some description here.
	"""
	
	def __init__(self, prop):
		self.BaseProblem = prop
		self.Psi = prop.psi
		self.PsiSize = prop.psi.GetData().size
		self.PotentialList = prop.Propagator.PotentialList
		if hasattr(self.BaseProblem.Propagator, "BasePropagator"):
			self.PotentialList = self.BaseProblem.Propagator.BasePropagator.PotentialList
	
		#Get time step and propagation time from baseproblem
		self.TimeStep = prop.TimeStep.real
		self.PropagationTime = prop.Duration
		self.PerturbControl = 0

		self.GetGoodBadRatio = False
		self.MaxYield = 0

		self.ControlVectorsMax = None
		self.ForwardSolutionMax = None


	def ApplyConfigSection(self, config):
		"""
		"""
	
		self.Config = config.Config
		self.ConfigSection = config

		#Get the list of control function potentials
		self.ControlFunctionNamesList = config.control_functions
		self.NumberOfControls = len(self.ControlFunctionNamesList)

		# Shortcuts to the salient control parameters
		self.ControlCutoff = config.control_cutoff
		self.MaxIterations = config.max_iterations
		self.YieldRequirement = config.yield_requirement

		if hasattr(config, "time_grid_size"):
			self.TimeGridSize = config.time_grid_size
			self.TimeGridResolution = self.PropagationTime / self.TimeGridSize
		else:
			self.TimeGridResolution = self.TimeStep
			self.TimeGridSize = int(round(self.PropagationTime / self.TimeStep))

		#Initial control
		self.InitialControl = config.__dict__.get("initial_control", InitialControl.Constant)
		print "    Initial control: %s" % self.InitialControl.Name

		#Set control cutoff from time step
		err = 1e-6
		self.ControlCutoff = err / self.TimeStep**2

		#Check for the debug option
		self.Debug = False
		if hasattr(config, "debug"):
			self.Debug = config.debug

		#Do backward updates?
		self.UpdateBackward = hasattr(config, "update_backwards") and config.update_backwards or False
		print "    Update backwards: %s" % self.UpdateBackward

		#Use diagonal penalty matrix or not
		#self.PenaltyMatrixIsDiagonal = config.__dict__.get("penalty_matrix_diagonal", True)

		#Apply perturbation to control update?
		self.PerturbControl = config.__dict__.get("perturb_control", 0)
		print "    Applying control perturbation: %e" % self.PerturbControl

		#Which penalty method are we using?
		self.PenaltyMethod = config.penalty_method

		#Depending on penalty type, different settings apply
		self.EnergyPenalty = config.energy_penalty
		if self.PenaltyMethod == PenaltyMethod.Energy:
			pass
		elif self.PenaltyMethod == PenaltyMethod.Projection:
			self.GoodSpacePenalty = config.good_space_penalty
			self.BadSpacePenalty = config.bad_space_penalty
			self.ControlSpaceProjectionMatrices = config.control_projection_matrix
			self.GoodMatrixPath = config.good_matrix_path
			self.BadMatrixPath = config.bad_matrix_path
		else:
			raise Exception("Unknown penalty method!")
		print "    Control penalty method: %s" % PenaltyMethod.Names[self.PenaltyMethod]

		#Get control update mask function
		self.MaskFunction = config.__dict__.get("mask_function", lambda t : 1)


	def Setup(self):
		self.SetupCommon()
		self.SetupPenaltyMatrix()
		self.SetupInitialControl()


	def SetupCommon(self):
		self.TimeGrid = linspace(0, self.PropagationTime, self.TimeGridSize)

		#Keep the calculated J's (shows convergence)
		self.J = []
		self.Yield = []
		self.GoodBadRatio = []

		#Two matrices to hold previous backward and forward propagation solutions (for each time step)
		self.ForwardSolution = zeros((self.PsiSize, self.TimeGridSize), dtype = complex)
		self.BackwardSolution = zeros((self.PsiSize, self.TimeGridSize), dtype = complex)

		#The target state
		self.SetupTargetState()

		#Create a copy of the wavefunction to calculate H|psi>
		self.TempPsi = self.BaseProblem.psi.CopyDeep()
		self.TempPsi2 = self.BaseProblem.psi.CopyDeep()

		#Make a list of the control functions
		self.ControlFunctionList = []
		for controlFunc in self.ControlFunctionNamesList:
			for potential in self.PotentialList:
				if potential.Name == controlFunc:
					self.ControlFunctionList += [potential]

		#Check that we actually got some control functions
		if len(self.ControlFunctionList) == 0:
			raise Exception("Did not find any control functions!")

		#A vector to store the optimized control functions
		self.ControlVectors = ones((self.NumberOfControls, self.TimeGridSize), dtype = double)
		#for a in range(self.NumberOfControls):
		#	self.ControlVectors[a,:] *= self.ControlFunctionList[a].ConfigSection.strength

	
	def Run(self):
		"""
		Run optimal control algorithm, until desired yield or MaxIterations is reached.
		"""

		#
		# Main control loop
		# Iterate until converged or MaxIterations is reached
		#
		self.currIter = 0
		self.__YieldNotReached = True
		while (self.currIter < self.MaxIterations) & self.__YieldNotReached:
			
			# Forward propagation
			self.ForwardPropagation()

			# Project on target state and check yield.
			targetProjection = self.ComputeTargetProjection()
			currentYield = abs(targetProjection)**2
			self.Yield.append(currentYield)

			#Check if we have a new max yield. If so, store control and forward solution
			if currentYield > self.MaxYield:
				self.MaxYield = currentYield
				self.ControlVectorsMax = self.ControlVectors.copy()
				self.ForwardSolutionMax = self.ForwardSolution.copy()

			#Get value of cost functional
			curJ = self.ComputeCostFunctional(currentYield)
			self.J.append(curJ)
			norm = self.BaseProblem.psi.GetNorm()
			print "Yield = %.15f, J = %.15f, norm = %.15f" % (currentYield, curJ, norm)

			#Get value of good/bad ratio if applicable
			if self.GetGoodBadRatio:
				if self.PenaltyMethod == PenaltyMethod.Projection:
					good = self.ProjectionMatrixExpectationValue(whichMatrix="Good")
					bad = self.ProjectionMatrixExpectationValue(whichMatrix="Bad")
					self.GoodBadRatio.append(sqrt(bad/good))

				if self.PenaltyMethod == PenaltyMethod.ProjectionSparse:
					good = self.ProjectionMatrixExpectationValue(whichMatrix="Good")
					bad = numpy.dot(self.ControlVectors[0], self.ControlVectors[0])
					self.GoodBadRatio.append(sqrt(bad/good))

			#Desired yield reached? -> Finished
			if( currentYield > self.YieldRequirement ):
				print "Yield reached!"
				self.__YieldNotReached = False
				continue

			# Backward propagation
			self.BackwardPropagation(targetProjection)

			self.currIter += 1


	def ComputeNewControlFunctions(self, timeGridIndex, t, direction):
		"""
		Updates the control function based on backward propagation solution:

		    E(t) = - 1 / mu * imag(<etta|X1 + X2 + O(h) |psi>)

		where E is the new control, mu is the energy penalty, |etta> is the
		backward solution and X the control matrix. The number of higher order
		terms depend on the algorithm at hand (Krotov, Zhu-Rabitz or Degani)
		"""

		#Check if we should skip backward update
		if (direction == Direction.Backward) and (not self.UpdateBackward):
			for a in range(self.NumberOfControls):
				self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, timeGridIndex]

		else:

			self.SetupVectorB(timeGridIndex, direction)
			self.SetupMatrixM(timeGridIndex, direction)
			newControls = linalg.solve(self.M, self.b)

			if self.PerturbControl > 0:	
				newControls[0] += [(rand() - 0.5) * self.PerturbControl]

			#Update controls
			self.UpdateControls(newControls, timeGridIndex)
	

	def ComputeCostFunctional(self, currentYield):
		penalty = 0
		for a in range(self.NumberOfControls):
			controlVec = self.ControlVectors[a,:]
			penalty += self.TimeGridResolution * self.PenaltyMatrix.ExpectationValue(controlVec)
		return currentYield - penalty


	def ComputeTargetProjection(self):
		"""
		Project solution on target state
		"""
		#Note: Argument is complex conjugated
		return self.BaseProblem.psi.InnerProduct(self.TargetState)


	def UpdateControls(self, newControls, timeGridIndex):
		"""
		Update controls from list newControls
		"""
		maskValue = self.MaskFunction(self.TimeGrid[timeGridIndex])
		for a in range(self.NumberOfControls):
			self.ControlVectors[a, timeGridIndex] = maskValue * newControls[a]
			self.ControlFunctionList[a].ConfigSection.strength = maskValue * newControls[a]


	def ForwardPropagation(self):
		"""
		Propagate forward from 0 to T.
		"""
		
		print "Forward propagation, pass ", self.currIter

		self.SetupStep(Direction.Forward)
		self.InitializeControls(Direction.Forward)
		
		#Initial step
		if self.currIter > 0:
			self.ComputeNewControlFunctions(0, 0.0, Direction.Forward)

		for idx, t in enumerate(self.BaseProblem.Advance(self.TimeGridSize)):
			self.ForwardSolution[:, idx] = self.BaseProblem.psi.GetData()[:]
			if (idx < self.TimeGridSize - 1) and (self.currIter > 0):
				self.ComputeNewControlFunctions(idx+1, t, Direction.Forward)
			else:
				for a in range(self.NumberOfControls):
					self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, idx]
	
	
	def BackwardPropagation(self, targetProjection):
		"""
		Propagate backwards from T to 0. Terminal condition is Psi(T) = <F|Psi_forward(T)> |F>,
		where |F> is the desired final state.
		"""

		print "Backward propagation, pass ", self.currIter

		self.SetupStep(Direction.Backward)
		#self.InitializeControls(Direction.Backward)

		# Set the initial state
		self.BaseProblem.psi.GetData()[:] = targetProjection * self.TargetState.GetData()[:]
		#self.BaseProblem.psi.GetData()[:] = self.ForwardSolution[:,-1]
		
		#Update controls
		self.ComputeNewControlFunctions(self.TimeGridSize - 1, self.PropagationTime, Direction.Backward)

		for idx, t in enumerate(self.BaseProblem.Advance(self.TimeGridSize)):
			self.BackwardSolution[:, self.TimeGridSize - 1 - idx] = self.BaseProblem.psi.GetData()[:]
			if idx < (self.TimeGridSize - 1):
				self.ComputeNewControlFunctions(self.TimeGridSize - idx - 2, t, Direction.Backward)

	
	def SetupStep(self, direction):
		"""
		Set up propagator for forward or backward run
		"""

		if direction == Direction.Forward:
			self.BaseProblem.RestartPropagation(complex(self.TimeStep), 0.0, self.PropagationTime)

		elif direction == Direction.Backward:
			self.BaseProblem.RestartPropagation(-complex(self.TimeStep), self.PropagationTime, self.PropagationTime)

		#Apply initial condition
		self.BaseProblem.SetupWavefunction()

	
	def SetupPenaltyMatrix(self):
		"""
		Set up penalty matrix. This is either a diagonal matrix with size
		equal to the number of control time intervals, or a full matrix
		consisting of projections onto a "good" and "bad" control space.
		In case of diagonal matrix, we only store the diagonal. Else we
		set up the full matrix using "good"/"bad" projection matrices
		read from a user-specified file.
		"""
		
		if self.PenaltyMethod == PenaltyMethod.Energy:
			#self.PenaltyMatrix = numpy.ones(self.TimeGridSize)
			#self.PenaltyMatrix[:] *= self.EnergyPenalty * self.TimeGridResolution
			#diagEntries = numpy.ones(self.TimeGridSize) * self.EnergyPenalty * self.TimeGridResolution
			self.PenaltyMatrix = PenaltyMatrix((self.TimeGridSize))
			diagEntries = numpy.ones(self.TimeGridSize, dtype=double) * self.EnergyPenalty
			self.PenaltyMatrix.SetDiagonal(diagEntries)

		elif self.PenaltyMethod == PenaltyMethod.Projection:
			penaltyMatrixMem = self.TimeGridSize**2 * 8 / (1024)**2
			print "    Creating penalty matrix of shape [%s,%s] (~%s MB)" \
				% (self.TimeGridSize, self.TimeGridSize, penaltyMatrixMem) 
			#self.PenaltyMatrix = numpy.zeros((self.TimeGridSize, self.TimeGridSize))
			self.PenaltyMatrix = PenaltyMatrix((self.TimeGridSize, self.TimeGridSize))
			self.ProjectionMatrixGood = PenaltyMatrix((self.TimeGridSize, self.TimeGridSize))
			self.ProjectionMatrixBad = PenaltyMatrix((self.TimeGridSize, self.TimeGridSize))

			#Set up diagonal energy penalty
			energyPenalty = numpy.ones(self.TimeGridSize, dtype=double)
			energyPenalty[:] *= self.EnergyPenalty
			self.PenaltyMatrix.SetDiagonal(energyPenalty)

			#Load projection matrices for "good" and "bad" control spaces
			h5file = tables.openFile(self.ControlSpaceProjectionMatrices, "r")
			try:
				self.ProjectionMatrixGood[:] = h5file.getNode(self.GoodMatrixPath)[:]
				self.ProjectionMatrixBad[:] = h5file.getNode(self.BadMatrixPath)[:]
				self.PenaltyMatrix[:] += self.GoodSpacePenalty * self.ProjectionMatrixGood[:]
				self.PenaltyMatrix[:] += self.BadSpacePenalty * self.ProjectionMatrixBad[:]
			finally:
				h5file.close()

		elif self.PenaltyMethod == PenaltyMethod.ProjectionSparse:
			#self.PenaltyMatrix = numpy.zeros((self.TimeGridSize, self.TimeGridSize))
			self.PenaltyMatrix = PenaltyMatrix(self.TimeGridSize)
			self.ProjectionMatrixGood = PenaltyMatrix((self.TimeGridSize, self.TimeGridSize))

			#Set up diagonal energy penalty
			energyPenalty = numpy.ones(self.TimeGridSize, dtype=double)
			energyPenalty[:] *= self.EnergyPenalty
			self.PenaltyMatrix.SetDiagonal(energyPenalty)

			#Load projection matrices for "good" and "bad" control spaces
			h5file = tables.openFile(self.ControlSpaceProjectionMatrices, "r")
			try:
				self.ProjectionMatrixGood = h5file.getNode(self.GoodMatrixPath).copy()
			finally:
				h5file.close()



	def SetupTargetState(self):
		"""
		Set up the desired target state
		"""
		self.TargetState = self.BaseProblem.psi.CopyDeep()
		self.TargetState.Clear()

		if self.BaseProblem.Config.FinalState.type == "vector":
			self.TargetState.GetData()[self.BaseProblem.Config.FinalState.states] \
				= self.BaseProblem.Config.FinalState.population[:]

		elif self.BaseProblem.Config.FinalState.type == "function":
			grid = self.Psi.GetRepresentation().GetLocalGrid(0)
			func = self.BaseProblem.Config.FinalState.grid_function
			self.TargetState.GetData()[:] = func(self.BaseProblem.Config.FinalState, grid)

		elif self.BaseProblem.Config.FinalState.type == InitialConditionType.File:
			h5file = tables.openFile(self.BaseProblem.Config.FinalState.filename, "r")
			try:
				self.TargetState.GetData()[:] = h5file.getNode(self.BaseProblem.Config.FinalState.dataset)[:]
			finally:
				h5file.close()

		else:
			raise Exception("Unknown target specification")
		
		#Normalize target state
		self.TargetState.Normalize()

	
	def SetupInitialControl(self):
		"""
		"""
		T = self.PropagationTime
		omega = 0.08
		f = self.InitialControl()
		for a in range(self.NumberOfControls):
			self.ControlVectors[a, :] = [f(t, T, self.ConfigSection) for t in linspace(0, T, self.TimeGridSize)]
			self.ControlVectors[a,:] *= self.ControlFunctionList[a].ConfigSection.strength

	
	def InitializeControls(self, direction):
		for a in range(self.NumberOfControls):
			self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, direction]


	def ProjectionMatrixExpectationValue(self, whichMatrix = "Good", whichControl = 0):
		if whichMatrix == "Good":
			return self.ProjectionMatrixGood.ExpectationValue(self.ControlVectors[whichControl])
		elif whichMatrix == "Bad":
			return self.ProjectionMatrixBad.ExpectationValue(self.ControlVectors[whichControl])
		else:
			print "Unknown projection matrix: %s" % whichMatrix


	def MultiplyCommutatorAB(self, A, B, psi, tmpPsi, outPsi):
		"""
		[A,B] * |psi> = AB * |psi> - BA * |psi>. Result is returned
		in outPsi.
		"""

		tmpPsi.Clear()
		outPsi.Clear()

		#Calculate -BA * |psi> store in outPsi
		A.MultiplyPotential(psi, tmpPsi, 0, 0)
		B.MultiplyPotential(tmpPsi, outPsi, 0, 0)
		outPsi.GetData()[:] *= -1.0

		#Calculate AB * |psi> store in outPsi
		tmpPsi.Clear()
		B.MultiplyPotential(psi, tmpPsi, 0, 0)
		A.MultiplyPotential(tmpPsi, outPsi, 0, 0)


	def MultiplyCommutatorAAB(self, A, B, psi, tmpPsi1, tmpPsi2, outPsi):
		"""
		Multiply the commutator [A,[A,B]]] on psi. Result is returned in outPsi.
		All buffers are destroyed.

		Straightforward expansion of the commutator would require 12 matrix-vector 
		multiplications. However, we write:

		    P1 = A |psi>
		    P2 = B |psi>

		giving

		    [A,[A,B]] = AAB - 2*ABA +  BAA = AA*P2 - 2*AB*P1 + BA*P1

		where we now must do only 8 matrix-vector multiplications.
		"""
		#Clear buffers
		tmpPsi1.Clear()
		tmpPsi2.Clear()
		outPsi.Clear()

		#Calculate P1
		A.MultiplyPotential(psi, tmpPsi1, 0, 0)

		#Calculate BA*P1 and store in outPsi
		A.MultiplyPotential(tmpPsi1, outPsi, 0, 0)
		B.MultiplyPotential(outPsi, tmpPsi2, 0, 0)
		outPsi.GetData()[:] = tmpPsi2.GetData()[:]

		#Calculate -2*AB*P1 and add to outPsi
		tmpPsi2.Clear()
		B.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		A.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] -= tmpPsi1.GetData()[:]

		#Calculate P2
		tmpPsi1.Clear()
		B.MultiplyPotential(psi, tmpPsi1, 0, 0)

		#Calculate AA*P2
		tmpPsi2.Clear()
		A.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		A.MultiplyPotential(tmpPsi2, tmpPsi1,0, 0)
		outPsi.GetData()[:] += tmpPsi1.GetData()[:]


	def MultiplyCommutatorABA(self, A, B, psi, tmpPsi1, tmpPsi2, outPsi):
		"""
		Multiply the commutator [A,[B,A]]] on psi. Result is returned in outPsi.
		All buffers are destroyed. Since [A,[B,A]] = -[A,[A,B]], we just call
		on MultiplyCommutatorAAB and multiply the result by -1
		"""
		self.MultiplyCommutatorAAB(A, B, psi, tmpPsi1, tmpPsi2, outPsi)
		outPsi.GetData()[:] *= -1


	def MultiplyCommutatorABC(self, A, B, C, psi, tmpPsi1, tmpPsi2, outPsi):
		"""
		Multiply the commutator [A,[B,C]] on psi. Result is returned in outPsi.

		[A,[B,C]] = [A,BC-CB] = ABC - ACB - BCA + CBA 
		"""

		def ClearBuffers():
			tmpPsi1.Clear()
			tmpPsi2.Clear()

		#Clear buffers
		ClearBuffers()
		outPsi.Clear()

		#Calculate CBA and store in outPsi
		A.MultiplyPotential(psi, tmpPsi1, 0, 0)
		B.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		C.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] = tmpPsi1.GetData()

		#Calculate -BCA and add to outPsi
		ClearBuffers()
		A.MultiplyPotential(psi, tmpPsi1, 0, 0)
		C.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		B.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] -= tmpPsi1.GetData()[:]

		#Calculate -ACB 
		ClearBuffers()
		B.MultiplyPotential(psi, tmpPsi1, 0, 0)
		C.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		A.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] -= tmpPsi1.GetData()[:]

		#Calculate ABC
		ClearBuffers()
		C.MultiplyPotential(psi, tmpPsi1, 0, 0)
		B.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		A.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] = tmpPsi1.GetData()[:]


class PenaltyMatrix:
	def __init__(self, size):
		if size.__class__ == int:
			self.Rank = 1
		else:
			self.Rank = len(size)
		if not (self.Rank == 1 or self.Rank==2):
			raise Exception("Penalty matrix must be of either rank 1 or two!")
		self.Size = size
		self.Matrix = numpy.zeros(size, dtype=double)

	def __setitem__(self, k, value):
		self.Matrix[k] = value
	
	def __getitem__(self, k):
		return self.Matrix[k]

	def GetRow(self, idx):
		return self.Matrix[idx,:]

	def GetDiagonal(self):
		if self.Rank == 1:
			return self.Matrix[:]
		else:
			return diag(self.Matrix)

	def GetDiagonalElement(self, idx):
		if self.Rank == 1:
			return self.Matrix[idx]
		else:
			return self.Matrix[(idx,idx)]

	def SetDiagonal(self, data):
		if self.Rank == 1:
			self.Matrix[:] = data
		else:
			self.Matrix[:] = diag(data)

	def Dot(self, a):
		if self.Rank == 1:
			return self.Matrix * a	
		else:
			return numpy.dot(self.Matrix, a)
			#return numpy.dot(a, self.Matrix)

	def DotRow(self, idx, vec):
		if self.Rank == 1:
			result = 0
		else:
			result = numpy.dot(self.Matrix[idx, :idx], vec[:idx])
			result += numpy.dot(self.Matrix[idx, idx+1:], vec[idx+1:])
			#result = numpy.dot(self.Matrix[idx, :], vec[:])
		return result

	def ExpectationValue(self, vec):
		return numpy.dot(vec, self.Dot(vec))


class PenaltyMatrixSparse:
	def __init__(self, smallSize, largeSize):
		self.Size = size
		self.Matrix = numpy.zeros((largeSize, smallSize), dtype=double)
		self.Diagonal = numpy.zeros(largeSize, dtype=double)

	def __setitem__(self, k, value):
		self.Matrix[k] = value
	
	def __getitem__(self, k):
		return self.Matrix[k]

	def GetRow(self, idx):
		return self.Matrix[idx,:]

	def GetDiagonal(self):
		return self.Diagonal[:]

	def GetDiagonalElement(self, idx):
		return self.Diagonal[idx]

	def SetDiagonal(self, data):
		if self.Rank == 1:
			self.Matrix[:] = data
		else:
			self.Matrix[:] = diag(data)

	def Dot(self, a):
		if self.Rank == 1:
			return self.Matrix * a	
		else:
			return numpy.dot(self.Matrix, a)
			#return numpy.dot(a, self.Matrix)

	def DotRow(self, idx, vec):
		if self.Rank == 1:
			result = 0
		else:
			result = numpy.dot(self.Matrix[idx, :idx], vec[:idx])
			result += numpy.dot(self.Matrix[idx, idx+1:], vec[idx+1:])
			#result = numpy.dot(self.Matrix[idx, :], vec[:])
		return result

	def ExpectationValue(self, vec):
		return numpy.dot(vec, self.Dot(vec))



class Direction:
	Forward = 0
	Backward = -1


class PenaltyMethod:
	Energy = 0
	Projection = 1
	Names = ["Energy", "Projection"]


class InitialControl:

	class Constant:
		Name = "Constant"	
		def __call__(self, t, T, conf):
			return 1

	class SineSquaredPulse:
		Name = "SineSquaredPulse"
		def __call__(self, t, T, conf):
			omega = conf.omega
			exponent = conf.envelope_exponent
			return sin(pi * t / T)**exponent * cos(omega * t)

	class SineSquaredRamp:
		"""
		Pulse is ramped on/off with sine squared envelope,
		with constant envelope in between.
		"""
		Name = "SineSquaredRamp"
		def __call__(self, t, T, conf):
			rampOnDuration = conf.ramp_on_duration
			rampOffDuration = conf.ramp_off_duration
			omega = conf.omega

			if t < rampOnDuration:
				envelope = sin(pi * t / (2 * rampOnDuration))**2
			elif t > (T - rampOffDuration):
				envelope = sin(pi * (T - t) / (2 * rampOffDuration))**2
			else:
				envelope = 1.0

			return envelope * cos(omega * t)
			
