package cs475.RBM2;
import java.util.Random;

import cs475.RBM2.*;

public class RBMEnergy {
	private RBMParameters _parameters;
	private int _num_samples;
	private int _m, _n; // for coding ease
	private double[] _ph; // holds p(h_j=1)
	private Random _rand = new Random(0);
	
	// TODO: Add the required data structures and methods.
	private double sigma(double z){
		return 1.0/(1.0 + Math.pow(Math.E, -z));
	}
	
	private double dotProd(double[] a, double[] b){
		double prod = 0;
		for(int i=0; i < a.length; i++)
			prod += a[i]*b[i];
		return prod;
	}
	
	// Wi_ returns i-th row of W
	private double[] Wi_(int i){
		double[] ret = new double[_n];
		for(int j=0; j<_n; j++)
			ret[j] = this._parameters.weight(i+1, j+1);
		return ret;
	}
	
	// W_j returns j-th col of W
	private double[] W_j(int j){
		double[] ret = new double[_m];
		for(int i=0; i<_m; i++)
			ret[i] = this._parameters.weight(i+1, j+1);
		return ret;
	}
	
	private double randVal(double p){
		// p = probability
		double r = this._rand.nextDouble();
		return r < p ? 1 : 0;
	}
	
	public RBMEnergy(RBMParameters parameters, int numSamples) {
		this._parameters = parameters;
		this._num_samples = numSamples;
		this._m = this._parameters.numVisibleNodes();
		this._n = this._parameters.numHiddenNodes();
		this._ph = new double[this._n];
		this.calculateAllMarignals();
	}
	
	public void calculateAllMarignals(){
		double[] x = new double[_m];
		double[] h = new double[_n];
		int[] nh1 = new int[_n]; // counts h==1
		
		for (int j=0; j<_n; j++)
			nh1[j] = 0;
		
		// initialize first x
		for (int i = 0; i<_m; i++)
			x[i] = (i+1)%2==0 ? 1 : 0;
		
		for (int t=0; t<this._num_samples; t++){
			// generate h, based on x
			for (int j=0; j<this._n; j++){
				double[] w_j = this.W_j(j);
				double xw_j = this.dotProd(x, w_j);
				double dj = this._parameters.hiddenBias(j+1);
				double phj1 = sigma(xw_j+dj);
				h[j] = this.randVal(phj1);
				if(h[j]==1)
					nh1[j]++;
			}
			
			// generate x, based on h
			for (int i=0; i<this._m; i++){
				double[] wi_ = this.Wi_(i);
				double hwi_ = this.dotProd(h, wi_);
				double bi = this._parameters.visibleBias(i+1);
				double pxi1 = sigma(hwi_+bi);
				x[i] = this.randVal(pxi1);
			}
		}
		
		for (int j=0; j<_n; j++)
			this._ph[j] = (nh1[j]+0.0) / this._num_samples;
	}
	
	public double computeMarginal(int j) {
		// TODO: Add code here
		return this._ph[j-1];
	}
}
