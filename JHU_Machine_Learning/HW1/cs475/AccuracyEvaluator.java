package cs475;

import java.util.List;

public class AccuracyEvaluator extends Evaluator {

	@Override
	public double evaluate(List<Instance> instances, Predictor predictor) {
		// TODO Auto-generated method stub
		int nTruePredictions = 0;
		int nFalsePredictions = 0;
		for (Instance instance : instances) {
			Label trueLabel = instance.getLabel();
			// cannot evaluate if true label is unknown.
			if (trueLabel == null)
				continue;
			
			// get predicted label and compare
			Label predictedLabel = predictor.predict(instance);
			if(trueLabel.toString().equals(predictedLabel.toString()))
				nTruePredictions ++;
			else
				nFalsePredictions++;
		}
		if(nTruePredictions + nFalsePredictions == 0)
			throw new IllegalArgumentException("Not enough known-labelled instances.");
		
		double accuracy = ((double)nTruePredictions)/(nTruePredictions + nFalsePredictions); 
		return accuracy;
	}

}
