package cs475.RBM;

public class RBMEnergy {
	private RBMParameters _parameters;
	private int _iters;
	private double _eta;
	
	// TODO: Add the required data structures and methods.

	public RBMEnergy(RBMParameters parameters, int iters, double eta) {
		this._parameters = parameters;
		this._iters = iters;
		this._eta = eta;
	}
	
	public void learning() {
		// TODO: Add code here
	}

	public void printParameters() {
		//NOTE: Do not modify this function
		for (int i=0; i<_parameters.numVisibleNodes(); i++)
			System.out.println("b_" + i + "=" + _parameters.getVisibleBias(i));
		for (int i=0; i<_parameters.numHiddenNodes(); i++)
			System.out.println("d_" + i + "=" + _parameters.getHiddenBias(i));
		for (int i=0; i<_parameters.numVisibleNodes(); i++)
			for (int j=0; j<_parameters.numHiddenNodes(); j++)
				System.out.println("W_" + i + "_" + j + "=" + _parameters.getWeight(i,j));
	}
}
