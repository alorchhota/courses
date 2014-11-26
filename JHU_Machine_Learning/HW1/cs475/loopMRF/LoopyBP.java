
package cs475.loopMRF;

import java.util.Hashtable;

public class LoopyBP {

	private LoopMRFPotentials potentials;
	private int iterations;
	// add whatever data structures needed
	Hashtable<String, Double> messages = new Hashtable<String, Double>();
	public static final String FtoXDir = "f->x";
	public static final String XtoFDir = "x->f";
	
	
	public LoopyBP(LoopMRFPotentials p, int iterations) {
		this.potentials = p;
		this.iterations = iterations;
		this.train();
	}

	private String getKey(String direction, int x, int f, int xi){
		String sep = "->";
		switch(direction){
		case XtoFDir:
			return "x" + x + sep + "f" + f + "(" + xi + ")";
		case FtoXDir:
			return "f" + f + sep + "x" + x + "(" + xi + ")";
		default:
			return "InvalidDirection";	
		}
	}
	
	private void setMessage(String direction, int x, int f, int xi, double msg){
		String key = this.getKey(direction, x, f, xi);
		this.messages.put(key, msg);
	}
	
	private double getMessage(String direction, int x, int f, int xi){
		String key = this.getKey(direction, x, f, xi);
		return this.messages.get(key);
	}
	
	private void setFtoXMessage(int f, int x, int xi, double msg){
		this.setMessage(FtoXDir, x, f, xi, msg);
	}
	
	private double getFtoXMessage(int f, int x, int xi){
		return this.getMessage(FtoXDir, x, f, xi);
	}
	
	private void setXtoFMessage(int x, int f, int xi, double msg){
		this.setMessage(XtoFDir, x, f, xi, msg);
	}
	
	private double getXtoFMessage(int x, int f, int xi){
		return this.getMessage(XtoFDir, x, f, xi);
	}
	
	private void train(){
		int f = 0;
		int x = 0;
		
		int k = this.potentials.numXValues();
		int n = this.potentials.loopLength();
		
		// initialize ux1->fn+1(x1) = 1 and ux1->f2n(x1)=1 for all values of x1
		for(int xi=1; xi<=k; xi++){
			this.setXtoFMessage(1, n+1, xi, 1);
			this.setXtoFMessage(1, 2*n, xi, 1);
		}
		
		for(int iter=1; iter<=this.iterations; iter++){
			int nothing = 0;
			// part (a)
			for(int i=1; i<=n; i++){
				// compute ufn+i->x(i+1)%n
				f = n+i;
				//x = (i+1)%n;
				x = i==n ? 1 : i+1;
				int prevx = i;
				for(int xi=1; xi<=k; xi++){
					double msg = 0;
					for(int prevxi=1; prevxi<=k ; prevxi++){
						msg += (this.getXtoFMessage(prevx, f, prevxi) * this.potentials.potential(f, prevxi, xi));
					}	
					this.setFtoXMessage(f, x, xi, msg);
				}
				
				// compute ux1+i%n->fn+1+i%n
				f = n+1+i%n;
				x = 1+i%n;
				int prevf = n+i;
				int topf = x;
				for(int xi=1; xi<=k; xi++){
					double msg = this.getFtoXMessage(prevf, x, xi) * this.potentials.potential(topf, xi);
					this.setXtoFMessage(x, f, xi, msg);
				}	
				
			}
			
			// part (b)
			for(int i=n; i>=1; i--){
				// compute ufn+i->xi
				f = n+i;
				x = i;
				int prevx = x==n ? 1 : x+1;
				for(int xi=1; xi<=k; xi++){
					double msg = 0;
					for(int prevxi=1; prevxi<=k ; prevxi++){
						msg += (this.getXtoFMessage(prevx, f, prevxi) * this.potentials.potential(f, xi, prevxi));
					}	
					this.setFtoXMessage(f, x, xi, msg);
				}
				
				// compute uxi->fn+(i-2)%n+1
				//f = n+(i-2)%n+1;
				f = i==1 ? 2*n : n+i-1;
				x = i;
				int prevf = n+i;
				int topf = i;
				for(int xi=1; xi<=k; xi++){
					double msg = this.getFtoXMessage(prevf, x, xi) * this.potentials.potential(topf, xi);
					this.setXtoFMessage(x, f, xi, msg);
				}
				
			}
		}
				
	}
	
	public double[] marginalProbability(int x) {
		// TODO
		int k = this.potentials.numXValues();
		int n = this.potentials.loopLength();
		
		double[] mp = new double[k+1];
		mp[0] = 0;

		int topf = x;
		int leftf = x==1? 2*n : n+x-1;
		int rightf = n+x;
		
		// calculate all message
		for(int xi=1; xi<=k; xi++){
			mp[xi] = this.potentials.potential(topf, xi) * this.getFtoXMessage(leftf, x, xi) * this.getFtoXMessage(rightf, x, xi) ; 
		}
		
		// convert to probability
		double total = 0;
		for(int xi=1; xi<=k; xi++)
			total += mp[xi];
		for(int xi=1; xi<=k; xi++)
			mp[xi] /= total;
			
		return mp;
	}

}

