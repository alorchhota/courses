package cs475.test;

import cs475.FeatureVector;


class FeatureVectorTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		FeatureVector fv = new FeatureVector();
		fv.add(1, 100);
		fv.add(18, 18.18);
		fv.add(12, -12.12);
		
		System.out.println(
				fv.get(0) + " " + 
				fv.get(1) + " " + 
				fv.get(10) + " " +
				fv.get(12));
	}

}
