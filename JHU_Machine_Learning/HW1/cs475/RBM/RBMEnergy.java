package cs475.RBM;

public class RBMEnergy {
	private RBMParameters _parameters;
	private int _iters;
	private double _eta;
	
	private int _m, _n, _T; // for coding ease
	
	// TODO: Add the required data structures and methods.
	private double sigma(double z){
		return 1.0/(1.0 + Math.pow(Math.E, -z));
	}
	
	private double dotProd(double[] a, double[] b){
		double prod = 0;
		for(int i=0; i < a.length; i++)
			prod += a[i]*b[i];
		return prod;
	}
	
	// Wi_ returns i-th row of W
	private double[] Wi_(int i){
		double[] ret = new double[_n];
		for(int j=0; j<_n; j++)
			ret[j] = this._parameters.getWeight(i, j);
		return ret;
	}
	
	// W_j returns j-th col of W
	private double[] W_j(int j){
		double[] ret = new double[_m];
		for(int i=0; i<_m; i++)
			ret[i] = this._parameters.getWeight(i, j);
		return ret;
	}
	
	// returns t-th example x(t) vector
	private double[] vx(int t){
		double[] ret = new double[_m];
		for(int i=0; i<_m; i++)
			ret[i] = this._parameters.getExample(t, i);
		return ret;
	}
	
	// returns b vector
	private double[] vb(){
		double[] ret = new double[_m];
		for(int i=0; i<_m; i++)
			ret[i] = this._parameters.getVisibleBias(i);
		return ret;
	}
	
	// returns d vector
	private double[] vd(){
		double[] ret = new double[_n];
		for(int i=0; i<_n; i++)
			ret[i] = this._parameters.getHiddenBias(i);
		return ret;
	}
	
	private double energy(double[] x, double[] h){
		// calculate xWh
		double[] xw = new double[_n];
		for(int j=0; j<_n; j++){
			xw[j] = dotProd(x, W_j(j));
		}
		double xwh = dotProd(xw, h);
		
		// calculate bx
		double bx = dotProd(this.vb(), x);
		// calculate dh
		double dh = dotProd(this.vd(), h);
		
		double e = -xwh-bx-dh;
		return e;
	}
	
	private double p(double[] x, double[] h, double z){
		double e = energy(x,h);
		double p = 1.0 / z * Math.pow(Math.E, -e);
		return p;
	}
	
	private double[] int2bits(int n, int nbits){
		double[] bits = new double[nbits];
	    for (int i = 0; i < nbits; i++) {
	        bits[i] = (n & (1 << i)) != 0 ? 1 : 0;
	    }
	    return bits;
	}
	
	
	private double Z(){
		double M = Math.pow(2, _m);
		double N = Math.pow(2, _n);
		double z = 0;
		for(int i=0; i<M; i++){
			double[] x = int2bits(i, _m); 
			for(int j=0; j<N; j++){
				double[] h = int2bits(j, _n);
				z += Math.pow(Math.E, -energy(x, h));
			}
		}
		return z;
	}
	
	public RBMEnergy(RBMParameters parameters, int iters, double eta) {
		this._parameters = parameters;
		this._iters = iters;
		this._eta = eta;
	}
	
	private void initParams(){
		// init W
		for(int i=0; i<_m; i++){
			for(int j=0; j<_n; j+=2)
				this._parameters.setWeight(i, j, 0);
			for(int j=1; j<_n; j+=2)
				this._parameters.setWeight(i, j, 1);
		}
		
		// init b
		for(int i=0; i<_m; i+=2)
			this._parameters.setVisibleBias(i, 0);
		for(int i=1; i<_m; i+=2)
			this._parameters.setVisibleBias(i, 1);
		
		// init d
		for(int j=0; j<_n; j+=2)
			this._parameters.setHiddenBias(j, 0);
		for(int j=1; j<_n; j+=2)
			this._parameters.setHiddenBias(j, 1);
		
	}
	
	public void learning() {
		// TODO: Add code here
		this._m = this._parameters.numVisibleNodes();
		this._n = this._parameters.numHiddenNodes();
		this._T = this._parameters.numExamples();
		
		// parameter initialization
		initParams();
		
		double[][] W_new = new double[_m][_n];
		double[] b_new = new double[_m];
		double[] d_new = new double[_n];
		
		double M = Math.pow(2, _m);
		double N = Math.pow(2, _n);
		
		double[][] p_xh = new double[(int)M][(int)N];
		
		//System.out.println("===== Iter: 0 =====");
		//this.printParameters();
		
		for (int iter=0; iter<this._iters; iter++){
			// calculate Z and probabilities once in an iteration
			double z = 0;
			for(int i=0; i<M; i++){
				double[] x = int2bits(i, _m); 
				for(int j=0; j<N; j++){
					double[] h = int2bits(j, _n);
					p_xh[i][j] = Math.pow(Math.E, -energy(x, h));
					z += p_xh[i][j]; 
				}
			}
			
			for(int i=0; i<M; i++){
				for(int j=0; j<N; j++){
					p_xh[i][j] = p_xh[i][j] / z;
				}
			}
			
			/*// print intermediate outputs
			System.out.println("===== Iter: " + (iter+1) + " =====");
			for(int ti=0; ti<_m; ti++)
				System.out.print("x[" + (ti+1) +"]  ");
			for(int ti=0; ti<_n; ti++)
				System.out.print("h[" + (ti+1) +"]  ");
			System.out.println("p");
			
			for(int xi=0; xi<M; xi++){
				double[] x = int2bits(xi, _m);
				for(int hi=0; hi<N; hi++){
					double[] h = int2bits(hi, _n);
					
					for(int ti=0; ti<_m; ti++)
						System.out.print(x[ti] + "   ");
					for(int ti=0; ti<_n; ti++)
						System.out.print(h[ti] + "   ");
					System.out.println(p(x,h,z));
				}
			}
			*/
			
			// calculate d/dWij
			for(int i=0; i<_m; i++){
				for(int j=0; j<_n; j++){
					double dw1 = 0;
					double[] w_j = W_j(j);
					double dj = this._parameters.getHiddenBias(j);
					for(int t=0; t<_T; t++){
						double[] xt = vx(t);
						dw1 += sigma(dotProd(xt, w_j) + dj) * xt[i];
					}
					
					double dw2 = 0;
					for(int xi=0; xi<M; xi++){
						double[] x = int2bits(xi, _m);
						if (x[i]!=1) continue;
						for(int hi=0; hi<N; hi++){
							double[] h = int2bits(hi, _n);
							if(h[j]!=1) continue;
							dw2 += p(x, h, z);
						}
					}
					dw2 *= _T;
					
					W_new[i][j] = dw1 - dw2;
				}
			}
			
			// calculate d/dbi
			for(int i=0; i<_m; i++){
				double db1 = 0;
				for(int t=0; t<_T; t++){
					db1 += this._parameters.getExample(t, i);
				}
				
				double db2 = 0;
				for(int xi=0; xi<M; xi++){
					double[] x = int2bits(xi, _m);
					if (x[i]!=1) continue;
					for(int hi=0; hi<N; hi++){
						double[] h = int2bits(hi, _n);
						db2 += p(x, h, z);
					}
				}
				db2 *= _T;
				
				b_new[i] = db1 - db2;
			}
			
			// calculate d/ddj
			for(int j=0; j<_n; j++){
				double dd1 = 0;
				double[] w_j = W_j(j);
				double dj = this._parameters.getHiddenBias(j);
				for(int t=0; t<_T; t++){
					double[] xt = vx(t);
					dd1 += sigma(dotProd(xt, w_j) + dj);
				}
				
				double dd2 = 0;
				for(int xi=0; xi<M; xi++){
					double[] x = int2bits(xi, _m);
					for(int hi=0; hi<N; hi++){
						double[] h = int2bits(hi, _n);
						if(h[j]!=1) continue;
						dd2 += p(x, h, z);
					}
				}
				dd2 *= _T;
				
				d_new[j] = dd1 - dd2;
			}
			
			
			// update Wij, bi, dj
			for(int i=0; i<_m; i++){
				for(int j=0; j<_n; j++){
					double w_old = this._parameters.getWeight(i, j);
					this._parameters.setWeight(i, j, w_old + this._eta * W_new[i][j]);
				}
			}
			
			for(int i=0; i<_m; i++){
				double b_old = this._parameters.getVisibleBias(i);
				this._parameters.setVisibleBias(i, b_old + this._eta * b_new[i]);
			}
			
			for(int j=0; j<_n; j++){
				double d_old = this._parameters.getHiddenBias(j);
				this._parameters.setHiddenBias(j, d_old + this._eta * d_new[j]);
			}
			
			/*
			// print values to debug
			System.out.println("===== Iter: " + (iter+1) + " =====");
			// probabilities
			for(int j=0; j<N; j++){
				for(int i=0; i<M; i++){
					//System.out.println(j + ", " + i + ", " + p_xh[i][j]);
				}
			}
			
			// gradients
			for(int i=0; i<_m; i++){
				for(int j=0; j<_n; j++){
					System.out.println( "Grad W_" + i + "_"+ j +": " + W_new[i][j]);
				}
			}
			
			for(int i=0; i<_m; i++){
				System.out.println( "Grad b_" + i + ": " + b_new[i]);
				
			}
			
			for(int j=0; j<_n; j++){
				System.out.println( "Grad d_" + j + ": " + d_new[j]);
				
			}
			
			this.printParameters();
			*/
			
			
			
		}
		
		//System.out.println("===== Final output =====");
		
	}

	public void printParameters() {
		//NOTE: Do not modify this function
		for (int i=0; i<_parameters.numVisibleNodes(); i++)
			System.out.println("b_" + i + "=" + _parameters.getVisibleBias(i));
		for (int i=0; i<_parameters.numHiddenNodes(); i++)
			System.out.println("d_" + i + "=" + _parameters.getHiddenBias(i));
		for (int i=0; i<_parameters.numVisibleNodes(); i++)
			for (int j=0; j<_parameters.numHiddenNodes(); j++)
				System.out.println("W_" + i + "_" + j + "=" + _parameters.getWeight(i,j));
	}
}
