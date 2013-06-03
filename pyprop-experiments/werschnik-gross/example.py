#Import system modules
import sys
import os
from pylab import *
import numpy
import tables

#Load pyprop
pypropPath = os.environ["PYPROPPATH"]
sys.path.insert(1, pypropPath)
import pyprop
pyprop = reload(pyprop)

execfile("%s/examples/vector/optimalcontrol/example.py" % pypropPath)

#Load the project module
from libpotential import *

def SetupConfig(**args):
	#Decide which config file to use
	configFile = "config.ini"
	if "config" in args:
		configFile = args["config"]

	#Load the config file
	conf = pyprop.Load(configFile)

	#Modify the config
	if "imtime" in args:
		imtime = args["imtime"]
		propSection = conf.Propagation
		dt = abs(propSection.timestep)
		renormalize = False
		if imtime:
			dt = -1.0j * dt
			renormalize = True

		propSection.timestep = dt
		propSection.renormalization = renormalize

	if "amplitude" in args:
		amplitude = args["amplitude"]
		conf.DynamicPotential.amplitude = amplitude
		
	return conf


def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop


def FindGroundstate(**args):
	args['imtime'] = True
	prop = SetupProblem(**args)
	
	for t in prop.Advance(10):
		E = prop.GetEnergy()
		if pyprop.ProcId == 0:
			print "t = %f, E = %f" % (t, E)

	E = prop.GetEnergy()
	if pyprop.ProcId == 0:
		print "Ground State Energy = %f" % E

	SaveWavefunction("groundstate.h5", "/wavefunction", prop)

	return prop


def FindEigenstates(**args):
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	numEigs = prop.Config.Arpack.krylov_eigenvalue_count
	outFile = prop.Config.Arpack.hdf5_filename

	#Find eigenvectors
	solver.Solve()

	if pyprop.ProcId == 0:
		h5file = tables.openFile(outFile, "w")
		try:
			#Save eigenvectors
			for i in range(numEigs):
				prop.psi.GetData()[:] = solver.GetEigenvector(i)
				prop.psi.Normalize()
				h5file.createArray("/", "eigenvector%03i" % i, prop.psi.GetData())

			#Save eigenvalues
			h5file.createArray("/", "eigenvalues", solver.GetEigenvalues()[:])

			#Store config object
			h5file.setNodeAttr("/", "configObject", prop.Config.cfgObj)
		finally:
			h5file.close()
	
	return solver


def CalculateEigenBasisMatrixElements(prop, potential, eigFileName, **args):
	"""
	From eigenvectors of the Werschnik-Gross 1D potential, 
	calculate matrix elements in the eigenfunction basis. Result is
	stored as a dense matrix in a HDF5 file.
	"""

	outFile = args.has_key("outFile") and args["outFile"] or "out/wg_matrix.h5"

	#Set up buffer
	tmpPsi = prop.psi.Copy()
	tmpPsi2 = prop.psi.Copy()

	eigFile = tables.openFile(eigFileName)
	
	numberOfEigs = eigFile.root.eigenvalues[:].size

	M = zeros((numberOfEigs, numberOfEigs), dtype=complex)
	for i in range(numberOfEigs):
		#Load eigenvector i
		prop.psi.GetData()[:] = eigFile.getNode("/eigenvector%03i" % i)
		for j in range(i, numberOfEigs):
			#Load eigenvector j
			tmpPsi.GetData()[:] = eigFile.getNode("/eigenvector%03i" % j)

			#Multiply potential:  X|j>
			tmpPsi2.Clear()
			potential.MultiplyPotential(tmpPsi, tmpPsi2, 0, 0)
			
			#Calculate <i|X|j>
			M[i,j] = prop.psi.InnerProduct(tmpPsi2)
			if i != j:
				M[j,i] = numpy.conj(M[i,j])

	eigFile.close()

	#Store matrix
	h5file = tables.openFile(outFile, "w")
	try:
		h5file.createArray("/", "wg_matrix_elements", M)
	finally:
		h5file.close()


def RunStabilityExperiment(**args):
	args["perturbControl"] = True
		
	controlProblem = Setup(**args)

	controlProblem.Run()
	v1 = controlProblem.ControlVectors[0,:].copy()

	controlProblem = Setup(**args)
	controlProblem.Run()
	v2 = controlProblem.ControlVectors[0,:].copy()

	semilogy(linspace(0, controlProblem.PropagationTime, controlProblem.TimeGridSize), abs(v1 - v2))

	return v1, v2
