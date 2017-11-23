import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class obj {

	public static void main(String[] args) throws FileNotFoundException {
		
		
		File file = new File("C:/Users/farshad.toosi/Documents/smaple.txt");
		
		Scanner in = new Scanner (file);
		
		
		ArrayList <String> lines = new ArrayList<String> ();
		while(in.hasNextLine()) {
			lines.add(in.nextLine());
		}
		
		int v = lines.get(0).split(" ").length;
		
		double [][] mat = new double[lines.size()][v]; 
		
		for(int i=0;i<lines.size();i++) {
			String [] l = lines.get(i).split(" ");
			for(int j=0;j<v;j++) {
				mat[i][j] = Double.parseDouble(l[j]);
			}
		}
		
		
		MatrixObjectification pp = new MatrixObjectification(mat);
		Synchronization.laplacian(pp.getSimilarityMatrix());
		
		
		
		Synchronization sc = new Synchronization(pp.getSimilarityMatrix(), pp.getSimilarityMatrixZeros() );
		
		double [][] eignes = sc.getClusters(0.99);
		
		/*for(int i=0;i<eignes.length;i++) {
			for(int j=0;j<eignes[i].length;j++) {
				System.out.print(eignes[i][j]+" ");
			}System.out.println("   :"+i);
		}*/
		
		
	}
	
	
}
