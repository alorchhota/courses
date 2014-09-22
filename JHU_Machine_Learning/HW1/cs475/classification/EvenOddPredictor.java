package cs475.classification;

import java.util.List;
import java.util.Set;

import cs475.ClassificationLabel;
import cs475.FeatureVector;
import cs475.Instance;
import cs475.Label;
import cs475.Predictor;

public class EvenOddPredictor extends Predictor{

	private static final long serialVersionUID = 1L;

	public EvenOddPredictor(){
	}
	
	@Override
	public void train(List<Instance> instances) {
		// nothing to train, prediction depends only on the test features.
	}

	@Override
	public Label predict(Instance instance) {
		// TODO Auto-generated method stub
		if(instance == null)
			return null;
		
		double evenSum = 0, oddSum = 0;
		
		FeatureVector fv = instance.getFeatureVector();
		Set<Integer> featureIndexes = fv.getFeatureIndexes();
		for(Integer index : featureIndexes){
			if(index % 2 == 0){
				evenSum += fv.get(index);
			}
			else{
				oddSum += fv.get(index);
			}
		}
		
		//System.out.println("EvenSum: " + evenSum + ", OddSum: " + oddSum );
		
		ClassificationLabel label = null;
		if(evenSum >= oddSum)
			label = new ClassificationLabel(1);
		else
			label = new ClassificationLabel(0);
		
		return label;
	}
}
