#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class QuantumDotPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	double Omega;
	double B;
	double Beta;

	double OmegaSqr;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("omega", Omega);
		config.Get("b", B);
		config.Get("beta", Beta);

		OmegaSqr = sqr(Omega);
	}

	/*
	 * Called once every timestep
	 */
	void CurTimeUpdated()
	{
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		double pot = sqr(OmegaSqr) / (64.0 * B) * sqr(sqr(x));
		pot += -OmegaSqr/4.0 * sqr(x) + Beta * x * sqr(x);
		return pot;
	}
};

template<int Rank>
class DipoleOperator : public PotentialBase<Rank>
{
public:

	void ApplyConfigSection(const ConfigSection &config)
	{
	}

	/*
	 * Called for every grid point at every time step.
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		return x;
	}
};
