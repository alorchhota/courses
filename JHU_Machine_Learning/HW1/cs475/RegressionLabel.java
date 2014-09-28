package cs475;

import java.io.Serializable;

public class RegressionLabel extends Label implements Serializable {

	private static final long serialVersionUID = 1L;
	
	private double _label = -1;

	public RegressionLabel(double label) {
		// TODO Auto-generated constructor stub
		this._label = label;
	}
	
	public double getValue(){
		return this._label;
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return Double.toString(this._label);
	}

}
