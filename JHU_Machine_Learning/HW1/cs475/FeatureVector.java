package cs475;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Set;

public class FeatureVector implements Serializable {

	private static final long serialVersionUID = 1L;
	
	private HashMap<Integer, Double> _features = null;
	
	public FeatureVector(){
		this._features = new HashMap<Integer, Double>();
	}
	
	public void add(int index, double value) {
		// TODO Auto-generated method stub
		if(value != 0)
			this._features.put(index, value);
	}
	
	public double get(int index) {
		// TODO Auto-generated method stub
		if(!this._features.containsKey(index)){
			return 0;
		}
		
		return this._features.get(index);
	}
	
	public Set<Integer> getFeatureIndexes(){
		return this._features.keySet();
	}
	
	

}
