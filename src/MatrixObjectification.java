import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import Jama.Matrix;

public class MatrixObjectification {

	double [][] matrix;
	boolean similarity = false;
	double [][] SquareMatrix;
	double [] Eigen;
	double [][] SimilarityWithZeros;
	
	
	public MatrixObjectification(double [][] mat) throws FileNotFoundException {
		this.matrix = mat;
		concept(mat);
		this.SimilarityWithZeros = MatrixPreparation(mat); 
	}
	
	private void concept(double [][] mat) throws FileNotFoundException {
		File file = new File("C:/Users/farshad.toosi/Documents/PCA.CXT");
		PrintWriter out = new PrintWriter(file);
		int p = mat.length;
		int v = mat[0].length;
		out.println("B");
		out.println("PCA");
		out.println(p);
		out.println(v);
		out.println();
		for(int i=0;i<p;i++)
			out.println("a"+(i+1));
		for(int i=0;i<v;i++)
			out.println("b"+(i+1));
		for(int i=0;i<p;i++) {
			for(int j=0;j<v;j++) {
				if(mat[i][j]==1)
					out.print("X");
				else
					out.print(".");
			}
			out.println();
		}
		
		out.close();
	}
	
	public double [][] MatrixPreparation(double [][] mat) {
		int m = mat.length;
		int n = mat[0].length;
		double [][] Nmat = new double [m][n*2];
		for(int i=0;i<n;i++) {
			double cont = 0;
			for(int j=0;j<m;j++) {
				if(mat[j][i] == 1) {
					cont++;
				}
			}
			if(cont == 0)
				cont = 0;
			else
				cont = 1.0-(cont/(double)m);
			
			for(int j=0;j<m;j++) {
				if(mat[j][i]==1)
					Nmat[j][i] = cont;
				else
					Nmat[j][i] = 0;
			}
		}
		/////// now for the complementory variables.
		for(int i=n;i<n*2;i++) {
			double cont = 0;
			for(int j=0;j<m;j++) {
				if(mat[j][i-n] == 0) {
					cont++;
				}
			}
			if(cont == 0)
				cont = 0;
			else
				cont = 1.0-(cont/(double)m);
			
			for(int j=0;j<m;j++) {
				if((1-mat[j][i-n])==1)
					Nmat[j][i] =  cont;
				else
					Nmat[j][i] = 0;
			}
		}
		TOSTRING(Nmat);
		return Nmat;
	}
	
	
	public double [][] getSimilarityMatrix() {
		
		Matrix A = new Matrix(matrix);
		Matrix B = A.transpose();
		Matrix C = A.times(B);
		
		
		SquareMatrix = C.getArray();
		similarity = true;
		return SquareMatrix;
		
	}
	
	public double [][] getSimilarityMatrixZeros() {
		Matrix A = new Matrix(SimilarityWithZeros);
		Matrix B = A.transpose();
		Matrix D = A.times(B);
		D = Mult(A, A.getRowDimension(), A.getColumnDimension(), D);
		
		SimilarityWithZeros = D.getArray();
		return SimilarityWithZeros;
	}
	

	
	
	private Matrix Mult (Matrix A, int p, int v, Matrix C) {
		System.out.println(v+"   VV  "+getSimilarityMatrix().length+"  "+getSimilarityMatrix()[0].length+"  "+p+"  ");
		for(int n=0;n<p;n++) {
			for(int m=n+1;m<p;m++) {
				double coh = 0;
				for(int j=0;j<v;j++) {
					int a1 = (j % (v/2));
					if(this.matrix[n][a1] != this.matrix[m][a1]) {
						coh -= (1/(double)v);
					}
				}				
				C.set(n, m, C.get(n, m) + coh);
				C.set(m, n, C.get(m, n) + coh);
			}
		}
		
		TOSTRING(C.getArray());
		return C;
	}
	
	
	public double [] getEigenValues() {
		
		if(similarity) {
			Matrix A = new Matrix (SquareMatrix);
			Eigen = A.eig().getRealEigenvalues();
		}
		else {
			getSimilarityMatrix();
			Matrix A = new Matrix (SquareMatrix);
			Eigen = A.eig().getRealEigenvalues();
		}
		return Eigen;
	}
	
	public static double [] getEigenValuesSep(double [][] Sq) {
		
		double [] Eig;
			Matrix A = new Matrix (Sq);
			Eig = A.eig().getRealEigenvalues();
		
		return Eig;
	}
	
	
	public static ArrayList<Set<Integer>> Groups (double [][] mat, boolean prnt, int v) {
		
		
		ArrayList<Set<Integer>> group = new ArrayList<Set<Integer>>();
		ArrayList<Integer> Mem = new ArrayList<Integer>();
		int g = 2;
		boolean playedOff = false;
		Set<Integer> Gmem = new HashSet();
		
		for(int i=0;i<v;i++) {
			for(int j=i+1;j<v;j++) {
				if(mat[i][j]==1) {
					Mem.add(i);
					Mem.add(j);
					mat[i][j] = g;
					mat[j][i] = g;
					playedOff = true;
					for(int p=0;p<Mem.size();p++) {
						for(int q=0;q<v;q++) {
							if(mat[q][Mem.get(p)]==1) {
								mat[q][Mem.get(p)] = g;
								mat[Mem.get(p)][q] = g;
								Mem.add(q);
								Gmem.add(q);
								playedOff = true;
							}
						}
					}
				}
				if(playedOff) {
					playedOff = false;
					group.add(Gmem);
					Gmem = new HashSet();
					g++;
				}
			}
		}
		for(int i=0;i<v;i++) {
			if(!Mem.contains(i)) {
				Gmem.add(i);
				group.add(Gmem);
				Gmem = new HashSet();
				mat[i][i] = g;
				g++;
			}
		}
		
		// degree of coherency
		
		if(prnt) {
			group.forEach(gr -> {
				gr.forEach(m -> {
					System.out.print(m+" ");
				});System.out.println();
			});
			
		}
		
		
		return group;
		
	}
	
	
	public static double degreeEva(ArrayList<Set<Integer>> group, double [][] mat) {
		
		double inside=0;
		double outside=0;
		
		for(int i=0;i<group.size();i++) {
			ArrayList <Integer> e1 = new ArrayList<Integer> ();
			e1.addAll(group.get(i));
			double in = 0;
			double ind =0;
			for(int n=0;n<e1.size();n++) {
				for(int m=n+1;m<e1.size();m++) {
					if(e1.get(n) != e1.get(m)) {
						double kk = mat[e1.get(n)][e1.get(m)];
						if(kk<0)
							kk=0;
						in += kk;
						ind++;
					}
				}
			}
			if(e1.size() == 1) 
				inside += 1;
			else
				inside += in/ind;
		}
		
		///////outsider 
		for(int i=0;i<group.size();i++) {
			ArrayList <Integer> e1 = new ArrayList<Integer> ();
			e1.addAll(group.get(i));
			for(int j=i+1;j<group.size();j++) {
				ArrayList <Integer> e2 = new ArrayList<Integer> ();
				e2.addAll(group.get(j));
			double out = 0;
			double ind =0;
			for(int n=0;n<e1.size();n++) {
				for(int m=0;m<e2.size();m++) {
					if(e1.get(n) != e2.get(m)) {
						double kk = mat[e1.get(n)][e2.get(m)];
						if(kk<0)
							kk=0;
						out += kk;
						ind++;
					}
				}
			}
			if(ind==0)
				outside += 1;
			else
				outside += out/ind;
			}
		}
		
		
		double degree = 0;
		if(outside!=0)
			degree = inside/(outside);
		else
			degree = inside+1;
		return degree;
	}
	
	public void TOSTRING(double [][] mat) {
		System.out.println();
		for(int i=0;i<mat.length;i++) {
			for(int j=0;j<mat[0].length;j++) {
				System.out.print(mat[i][j]+" ");
			}System.out.println();
		}
	}
	
	
}
