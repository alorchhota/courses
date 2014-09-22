package cs475.classification;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import cs475.ClassificationLabel;
import cs475.Instance;
import cs475.Label;
import cs475.Predictor;

public class MajorityPredictor extends Predictor{

	private static final long serialVersionUID = 1L;

	private ClassificationLabel _majorityLabel = null;
	
	public MajorityPredictor(){
		this._majorityLabel = null;
	}
	
	@Override
	public void train(List<Instance> instances) {
		// TODO Auto-generated method stub
		
		HashMap<String, Integer> labelCounts = new HashMap<String, Integer>();
		int maxCount = 0;
		ArrayList<String> maxLabels = new ArrayList<String>(); 
		
		for(Instance instance : instances){
			Label l = instance.getLabel();
			if(l == null)
				continue;
			
			String label = l.toString();
			int labelCount = labelCounts.containsKey(label) ? labelCounts.get(label)+1 : 1;
			labelCounts.put(label, labelCount);
			
			if(labelCount > maxCount){
				maxLabels.clear();
				maxLabels.add(label);
				maxCount = labelCount;
			}
			else if (labelCount == maxCount && !maxLabels.contains(label)){
				maxLabels.add(label);
			}
		}
		
		// randomly pick a max label
		Random rand = new Random();
		int selectedIndex = rand.nextInt(maxLabels.size());
		int majorityLabel = Integer.parseInt(maxLabels.get(selectedIndex));
		this._majorityLabel = new ClassificationLabel(majorityLabel);
		
	}

	@Override
	public Label predict(Instance instance) {
		// TODO Auto-generated method stub
		if(instance == null)
			return null;
		
		return this._majorityLabel;
	}
	

}
