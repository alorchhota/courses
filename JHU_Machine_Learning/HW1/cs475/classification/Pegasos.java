package cs475.classification;

import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;

import cs475.ClassificationLabel;
import cs475.FeatureVector;
import cs475.Instance;
import cs475.Label;
import cs475.Predictor;

public class Pegasos extends Predictor{

	private static final long serialVersionUID = 1L;
	
	private int D;				// Dimension
	private int _iter;			// #iterations
	private double _lambda;	// n0 (eta_0) for AdaGrad
	private double[] _w;
	
	
	DecimalFormat _fmt = new DecimalFormat("#.####");
	
	private double ll(boolean truthValue){
		return truthValue ? 1 : 0;
	}
	
	private double dotProd(double[] w, double[] x){
		double dp = 0;
		for(int i=0; i<w.length; i++){
			dp += w[i]*x[i];
		}
		return dp;
	}
	
	private int getLabel(Instance inst){
		int l = ((ClassificationLabel)inst.getLabel()).getValue();
		return l==0 ? -1 : l;
	}
	
	private double[] getFeatures(Instance inst){
		double[] x = new double[D];
		FeatureVector fv = inst.getFeatureVector();
		Set<Integer> fIndexes = fv.getFeatureIndexes();
		for(int fi : fIndexes){
			if(fi <= this.D)
				x[fi-1] = fv.get(fi);	// feature index is 1-based
		}
		return x;
	}
	
	public Pegasos(int iter, double lambda){
		this._iter = iter;
		this._lambda = lambda;
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
		for(int i=0; i<this.D; i++){
			w[i] = 0;
		}
		
		// populate some variables
		int N = instances.size();
		
		
		
		// pegasos
		for(int iter=0; iter < this._iter * N; iter++){
			int t = iter + 1;
			Instance curInstance = instances.get(iter%N);
			
			// extract y from sample
			int y = this.getLabel(curInstance);
			
			// extract x from sample
			double[] x = this.getFeatures(curInstance);
			
			
			double multiplier1 = 1.0 - 1.0/t; 
			double wx = dotProd(w, x);
			double llval = ll(y*wx < 1);
			double multiplier2 = 1.0 / (this._lambda * t) * llval * y;
			
			for(int d=0; d<this.D; d++){
				// update w
				w[d] =  multiplier1 * w[d] + multiplier2 * x[d];
			}
			
		}
		
		// save weights of the model
		this._w = new double[this.D];
		for(int d=0; d<this.D; d++)
			this._w[d] = w[d];
		
//		for(int d=0; d<this.D; d++)
//			System.out.print(this._fmt.format(w[d]) + "\t");
//		System.out.println();
	}

	@Override
	public Label predict(Instance instance) {
		if(instance == null)
			return null;
		
		// extract x from sample
		double[] x = this.getFeatures(instance);
		
		double wx = this.dotProd(this._w, x);
		int l =  wx>=0 ? 1 : 0;
		
		//System.out.println(l + "\t" + wx);
		
		return new ClassificationLabel(l);
	}
}
