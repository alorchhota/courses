package cs475.classification;

import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;

import cs475.ClassificationLabel;
import cs475.FeatureVector;
import cs475.Instance;
import cs475.Label;
import cs475.Predictor;

public class LogisticRegressionSGD extends Predictor{

	private static final long serialVersionUID = 1L;
	
	private int D;				// Dimension
	private int _iter;			// #iterations
	private double _eta0;	// n0 (eta_0) for AdaGrad
	private double[] _w;
	
	
	DecimalFormat _fmt = new DecimalFormat("#.####");
	
	private double g(double z){
		return 1.0 / (1+ Math.pow(Math.E, -z));
	}
	
	public LogisticRegressionSGD(int iter, double eta0){
		this._iter = iter;
		this._eta0 = eta0;
		this._w = null; 
	}
	
	private int getFeautreDimension(List<Instance> instances){
		// highest feature index is the total number of Features
		int maxFeatureIndex = 0;
		for(Instance ins : instances){
			Set<Integer> featureIndexes = ins.getFeatureVector().getFeatureIndexes();
			for(Integer fi : featureIndexes){
				if(fi > maxFeatureIndex){
					maxFeatureIndex = fi;
				}
			}
		}
		
		return maxFeatureIndex;
	}
		
	@Override
	public void train(List<Instance> instances) {
		this.D = this.getFeautreDimension(instances);
		
		// initialize w and eta
		double[] w = new double[this.D];
		double[] eta = new double[this.D];
		double[] I = new double[this.D];
		for(int i=0; i<this.D; i++){
			w[i] = 0;
			eta[i] = this._eta0;
			I[i] = 1;
		}
		
		// populate sum variables
		int N = instances.size();
		
		
		
		// stochastic gradient descent iteration
		double[] sum_f2 = new double[this.D];	// sum of f^2
		for(int i=0; i<this.D; i++){
			sum_f2[i] = 0;
		}
		
		for(int iter=0; iter < this._iter; iter++){
			Instance curInstance = instances.get(iter%N);
			
			// extract y from sample
			int y = ((ClassificationLabel)curInstance.getLabel()).getValue();
			
			// extract x from sample
			double[] x = new double[D];
			FeatureVector fv = curInstance.getFeatureVector();
			Set<Integer> fIndexes = fv.getFeatureIndexes();
			for(int fi : fIndexes){
				x[fi-1] = fv.get(fi);	// feature index is 1-based
			}
			
			// calculate gradient of likelihood for each dimension
			double[] dl = new double[this.D];
			double wx = 0;
			for(int d=0; d<this.D; d++){
				wx += w[d]*x[d];
			}
			for(int d=0; d<this.D; d++){
				dl[d] = y * g(-wx) * x[d] + (1-y) * g(wx) * (-x[d]);
			}
			
			// update sum_f2
			for(int d=0; d<this.D; d++){
				sum_f2[d] = sum_f2[d] + dl[d]*dl[d];
			}
			
			// update eta
			for(int d=0; d<this.D; d++){
				eta[d] = this._eta0 / (Math.sqrt(I[d] + sum_f2[d]));
			}
			
			// update w
			for(int d=0; d<this.D; d++){
				w[d] = w[d] + eta[d] * dl[d];
			}
			
			// print w
			System.out.print("Iteration#" + (iter+1) + ":\t");
			for(int d=0; d<this.D; d++){
				System.out.print(this._fmt.format(w[d]) + "\t");
			}
			System.out.println();
		}
		
		// save weights of the model
		this._w = new double[this.D];
		for(int d=0; d<this.D; d++)
			this._w[d] = w[d];
		
		for(int d=0; d<this.D; d++)
			System.out.print(this._fmt.format(w[d]) + "\t");
		System.out.println();
	}

	@Override
	public Label predict(Instance instance) {
		if(instance == null)
			return null;
		
		// extract x from sample
		double[] x = new double[D];
		FeatureVector fv = instance.getFeatureVector();
		Set<Integer> fIndexes = fv.getFeatureIndexes();
		for(int fi : fIndexes){
			x[fi-1] = fv.get(fi);	// feature index is 1-based
		}
		
		// y = g(wx)
		double wx = 0;
		for(int d=0; d<D; d++){
			wx += this._w[d] * x[d];
		}
		
		double y = g(wx);
		int l =  y>=0.5 ? 1 : 0;
		
		System.out.println(l + "\t" + wx);
		
		return new ClassificationLabel(l);
	}
}
