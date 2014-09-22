package cs475;

import java.io.Serializable;

public class ClassificationLabel extends Label implements Serializable {

	private static final long serialVersionUID = 1L;
	
	private int _label = -1;
	
	public ClassificationLabel(int label) {
		// TODO Auto-generated constructor stub
		this._label = label;
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return Integer.toString(this._label);
	}

}
