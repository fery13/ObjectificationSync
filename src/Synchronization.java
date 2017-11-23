import java.util.ArrayList;
import java.util.Set;

import Jama.Matrix;

public class Synchronization {

	double [][] similarity;
	double [][] similarityZeros;
	//double [] eigen;
	int Main_Nu_groups;
	double strength;
	public Synchronization (double [][] sim, double [][] simZero) {
		this.similarity = sim;
		this.similarityZeros = simZero;
		for(int i=0;i<sim.length;i++) {
			for(int j=i;j<sim.length;j++) {
				if(sim[i][j]!=0)
					strength += sim[i][j];
			}
		}
	}
	
	
	public double [][]  getClusters (double tr) {
		ArrayList<int[]> Clusters = new ArrayList<int[]> ();
		double round = 300;
		int time = similarity.length*150;
		int v = similarity.length;
		double [][] no_clusters = new double [time][v]; 
		
		double [][][] sync = new double [time][v][v];
		
		//*************** The original matrix
		{
			Matrix G = new Matrix(laplacian (this.similarity));
			System.out.println(toString());
			//ArrayList<Set<Integer>> group = MatrixObjectification.Groups(this.similarity, true);
			double [] temp = G.eig().getRealEigenvalues();
			int count =0;
			for(int i=0;i<temp.length;i++) {
				if(Math.abs(temp[i])<0.01)
					count++;
			}
			Main_Nu_groups = count;
			if(Main_Nu_groups  == 0)
				Main_Nu_groups  = 1;
			System.out.println(count+"  Number of groups in original");
		}
		//***********************************
		
		
		
		for(int r=0; r<round; r++) {
		
				double [][] polars = new double [time][v];
				for(int i=0;i<v;i++) {
					polars[0][i] = Math.random()*Math.PI*2;
				}
				
				for(int i=0;i<v;i++) {
					for(int j=0;j<v;j++) {
						sync[0][i][j] += (Math.PI<Math.abs(polars[0][i]- polars[0][j])) ? (2*Math.PI- Math.abs(polars[0][i]- polars[0][j])) : Math.abs(polars[0][i]- polars[0][j]);
					}
				}
				/// Start Time
				for(int t=1; t<time; t++) {
					for(int i=0;i<v;i++) {
						for(int j=0;j<v;j++) {
							double coupling = similarityZeros[i][j];
							if (coupling<0) 
								coupling =0;
							polars[t][i] +=  sync(polars[t-1][i], polars[t-1][j], coupling, (double) (1.0/(strength)));
							polars[t][i] = within_a_circle(polars[t][i]);
						}
						polars[t][i] += polars[t-1][i];
						polars[t][i] = within_a_circle(polars[t][i]);						
					}
					
					//Similarities
					for(int i=0;i<v;i++) {
						for(int j=0;j<v;j++) {
							sync[t][i][j] += (Math.PI<Math.abs(polars[t][i]- polars[t][j])) ? (2*Math.PI-Math.abs(polars[t][i]- polars[t][j])) : Math.abs(polars[t][i]- polars[t][j]);
						}
					}
					////////////////
				}
				
		}
		
		for(int t=0; t<time; t++) {
			for(int i=0;i<v;i++) {
				for(int j=0;j<v;j++) {
					sync[t][i][j] /= round;
					sync[t][i][j] = (tr>Math.cos(sync[t][i][j]) ? 0 :1);
				}
			}
			
			Matrix G = new Matrix(laplacian (sync[t]));
			double [] temp = G.eig().getRealEigenvalues();
			no_clusters[t] = temp;
			int count =0;
			for(int i=0;i<no_clusters[t].length;i++) {
				if(Math.abs(no_clusters[t][i])<0.001)
					count++;
			}
			
			
				
			ArrayList<Set<Integer>> group = MatrixObjectification.Groups(sync[t], true,v);
			MatrixObjectification.degreeEva(group, similarity);
			//System.out.println("------------   "+count+"   ------------------   "+MatrixObjectification.degreeEva(group, sync[t]));
			System.out.println("# Objects:	"+count+"	"+MatrixObjectification.degreeEva(group, similarity));
			
			
			if(Main_Nu_groups == count) {
				break;
			}
		}	
		
		return no_clusters;
	}
	
	
	
	private static double sync(double A, double B, double link, double coupling) {
		double synced = 0;
		synced = link * Math.sin((B-A)) * coupling;
		return synced;
	}
	
	
	private static double within_a_circle (double ang) {
		return ((Math.abs(ang)>Math.PI*2) ?  (ang % (Math.PI*2)) : (ang));
	}
	
	public static double [][] laplacian (double [][] mat) {
		double [][] lap = new double [mat.length][mat[0].length];
		
		for(int i=0;i<mat.length;i++) {
			int deg =0;
			for(int j=0;j<mat[0].length;j++) {
				if(mat[i][j]!=0 && i!=j) {
					lap[i][j] = mat[i][j]*(-1);
					deg += mat[i][j];
				}
			}
			lap[i][i]=deg;
		}
		
		return lap;
	}
	
	public String toString() {
		String temp = "";
		for(int i=0;i<this.similarity.length;i++) {
			for(int j=0;j<this.similarity.length;j++) {
				temp = temp +"  "+(int)this.similarity[i][j];
				
			}temp = temp +"\n";
		}
		
		return temp;
	}
	
	public static String ToString(double [][] mat) {
		String temp = "";
		for(int i=0;i<mat.length;i++) {
			for(int j=0;j<mat.length;j++) {
				temp = temp +" "+(int)mat[i][j];
				
			}temp = temp +"\n";
		}
		
		return temp;
	}
	
	
	
}
