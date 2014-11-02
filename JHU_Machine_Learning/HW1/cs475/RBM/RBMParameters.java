package cs475.RBM;

import java.io.BufferedReader;
import java.io.FileReader;
//NOTE: Do not modify this class
public class RBMParameters {

	private double[][] _weights;
	private double[][] _examples;
	private double[] _visibleBias;
	private double[] _hiddenBias;
	private int m, n, T;

	public RBMParameters(String data_file) {
		int line_num = 1;
		try {
			BufferedReader br = new BufferedReader(new FileReader(data_file));
			String line = br.readLine();
			String[] ar = line.split("\\s+");
			m = Integer.parseInt(ar[0]);
			n = Integer.parseInt(ar[1]);
			T = Integer.parseInt(ar[2]);
			_weights = new double[m][n];
			_examples = new double[T][m];
			_visibleBias = new double[m];
			_hiddenBias = new double[n];
			
			for(int i=0; i<m; i++) { 
				_visibleBias[i] = 0.0;
			}
			for(int i=0; i<n; i++) { 
				_hiddenBias[i] = 0.0;
			}
			for(int i=0; i<m; i++) { 
				for(int j=0; j<n; j++) { 
					_weights[i][j] = 0.0;
				}
				for(int t=0; t<T; t++) { 
					_examples[t][i] = 0.0;
				}
			}

			for (int t=0; t<T; t++) {
				line_num++;
				line = br.readLine();
				ar = line.split("\\s+");
				//System.out.println(ar[0]);
				if(ar.length == m) {
					for(int i=0; i<m; i++) 
						_examples[t][i] = Double.parseDouble(ar[i]);
				}
				else {
					continue;
				}
			}
			br.close();
		} catch(Exception e) {
			throw new RuntimeException("error while reading parameter file: " + data_file
				+ " on line " + line_num + " [" + e.getMessage() + "]");
		}
	}
	
	public int numHiddenNodes() { return n; }

	public int numVisibleNodes() { return m; }

	public int numExamples() { return T; }
	
	public double getExample(int t, int index) {
		if(t < 0 || t >= T)
			throw new RuntimeException("given t=" + t + ", illegal value for example: " + t);
		if(index < 0 || index >= m)
			throw new RuntimeException("given m=" + m + ", illegal value for index: " + index);
		return _examples[t][index];
	}

	public double getVisibleBias(int index) {
		if(index < 0 || index >= m)
			throw new RuntimeException("given m=" + m + ", illegal value for index: " + index);
		return _visibleBias[index];
	}
	
	public double getHiddenBias(int index) {
		if(index < 0 || index >= n)
			throw new RuntimeException("given n=" + n + ", illegal value for index: " + index);
		return _hiddenBias[index];
	}
	
	public double getWeight(int i, int j) {
		if(i < 0 || i >= m)
			throw new RuntimeException("given m=" + m + ", illegal value for i: " + i);
		if(j < 0 || j >= n)
			throw new RuntimeException("given n=" + n + ", illegal value for j: " + j);
		return _weights[i][j];
	}

	public boolean setVisibleBias(int index, double val) {
		if(index < 0 || index >= m) {
			System.out.println("Invalid index");
			return false;
		}
		_visibleBias[index] = val;
		return true;
	}
	
	public boolean setHiddenBias(int index, double val) {
		if(index < 0 || index >= n) {
			System.out.println("Invalid index");
			return false;
		}
		_hiddenBias[index] = val;
		return true;
	}
	
	public boolean setWeight(int i, int j, double val) {
		if(i < 0 || i >= m) {
			System.out.println("Invalid i");
			return false;
		}
		if(j < 0 || j >= n) {
			System.out.println("Invalid j");
			return false;
		}
		_weights[i][j]=val;
		return true;
	}

	
}
