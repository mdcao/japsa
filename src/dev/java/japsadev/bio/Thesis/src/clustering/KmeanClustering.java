package clustering;

import java.util.ArrayList;
import clustering.PairDistance.EditDistanceResult;
import clustering.GettingTreadsFromFasta;

/**
 * @author buvan.suji
 *
 */

public class KmeanClustering {


	/**
	 * @param args
	 */

	static int count1,count2,count3;
	static int d[];
	static int k[][];
	static int tempk[][];
	static int m[];
	static double diff[];	
	static int n;	
	final static int p =2;//number of clusters
	static double max = 0;
	static double temp1 = 0;
	static int index1;
	static int index2;	
	static ArrayList<Integer> t = new ArrayList<Integer>() ;
	
	
	static ArrayList<String> descript = new ArrayList<String>();
	static ArrayList<String> reads = new ArrayList<String>();
	static int Nreads;
	static String FileName;
	static int NumberElements;
	static int ClustElements;
	static ArrayList<String> distList1 = new ArrayList<String>();
	static ArrayList<ArrayList<String>> distList2 = new ArrayList<ArrayList<String>>();
	
	
	public static void Clustering() throws Exception{
		GettingTreadsFromFasta.DestReads();
		descript = GettingTreadsFromFasta.GetRname();
		reads = GettingTreadsFromFasta.GetTReads();
		Nreads = GettingTreadsFromFasta.NumberReads();
		FileName = GettingTreadsFromFasta.GetFileName();
		NumberElements = GettingTreadsFromFasta.SeqLength();
		ClustElements = GettingTreadsFromFasta.NelementsClustering();			
		
		double[][] table = new double[Nreads][Nreads];	
		
		for (int i=0; i<Nreads; ++i){
			for(int j=0;j<Nreads;++j){
				table[i][j]=0;
			}
		}

		for (int i = 0; i < reads.size(); i++) {
			for (int j = i + 1; j < reads.size(); j++) {				
				String x = reads.get(i);
				String y = reads.get(j);				
				EditDistanceResult result = PairDistance.compute(x, y);
				double normdist = ((double) result.getDistance())
						/ Math.max(x.length(), y.length());				
				table[i][j]=normdist;				
			}			
		}
		
		n=Nreads;
		d=new int[n];
		
		 //name of the reads		
		for(int i=0;i<n;++i){
			d[i]=i;
		}	
		
		 //Initialising arrays 
		k=new int[p][n];
		tempk=new int[p][n];		
		m=new int [p];
		diff=new double[p];		
		
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				temp1 = table[i][j];
				if(max < temp1){
					max = temp1;
					index1=i;
					index2=j;
				}
			}
		}			
		
		 //Initializing m 
		m[0] = index1;
		m[1] = index2;		
		
		int temp=0;
		int flag=0;
		
		do{
			for(int i=0;i<p;++i){
				for(int j=0;j<n;++j){
					k[i][j]=-1;					
				}
			}
			
			for(int i=0;i<n;++i){
				temp=NewCluster(d[i], table);					
				if(temp==0){
					k[temp][count1++]=d[i];					
				}
				else if(temp==1){
					k[temp][count2++]=d[i];					
				}					 
			}
			
			NewMean(table); // call to method which will calculate mean at this step.
			flag = VerifyEqual(); // check if terminating condition is satisfied.
			if(flag!=1){
				/*Take backup of k in tempk so that you can check for equivalence in next step*/
				for(int i=0;i<p;++i){
					for(int j=0;j<n;++j){
						tempk[i][j]=k[i][j];
					}
				}				
			}
			
			System.out.println("\nClusters Results:");
			for(int i=0;i<p;++i){
				System.out.print("C"+(i+1)+"{ ");
				for(int j=0;k[i][j]!=-1 && j<n-1;++j)
				System.out.print(k[i][j]+" ");
				System.out.println("}");
			}
			//end of for loop
			
			System.out.println("\nCentroid of the clusters are: ");
			for(int i=0;i<p;++i){
				System.out.print("m"+(i+1)+"="+m[i]+"  ");
			}	
			count1=0;count2=0;count3=0;			
		}
		while(flag==0);
		
		System.out.println("\n\n\nThe Final Clusters By Kmeans are as follows: ");
		for(int i=0;i<p;++i){
			System.out.print("C"+(i+1)+"{ ");
			for(int j=0;k[i][j]!=-1 && j<n-1;++j){
				//System.out.print(k[i][j]+" ");
				System.out.println(descript.get(k[i][j]));				
			}
			
			System.out.println("}");
		}		
	}
	

	
	public static int NewCluster(int a, double list[][]){ 		
		//int temp1=0;		
		for(int i=0;i<p;++i){
			if(a<m[i]){
				diff[i] = list[a][m[i]];
			}
			else if(a>m[i]){
				diff[i] = list[m[i]][a];
			}
			else{
				diff[i]=0;
			}			
		}
		
		int val=0;
		double temp=diff[0];
		for(int i=0;i<p;++i){
			if(diff[i]<temp){
				temp=diff[i];
				val=i;
			}
			
		}//end of for loop
		
		return val;
	}
	
	static void NewMean(double list [][]){
		
		for(int i=0;i<p;++i){
			m[i]=0; // initializing means to 0			
		}
		
		//int cnt=0;
		for(int i=0;i<p;++i){			
			
			for(int j=0;j<n-1;++j){								
				if(k[i][j]!=-1){
					t.add(k[i][j]);
				}
			}
				
			ArrayList<Double> s = new ArrayList<Double>();				
				
			for(int x1 = 0; x1<t.size();++x1){
				double sum = 0;					
				for(int x2=0;x2<t.size();++x2){
					if(t.get(x1)<t.get(x2)){
						sum = sum+(list[t.get(x1)][t.get(x2)]*list[t.get(x1)][t.get(x2)]);
					}
					else if(t.get(x1)>t.get(x2)){
						sum = sum+(list[t.get(x2)][t.get(x1)]*list[t.get(x2)][t.get(x1)]);
					}
					else{
						sum = sum+0;
					}	
				}
					//System.out.println(sum);
					
					s.add(sum/t.size());					
			}
			//System.out.println(s.size());
			double min=s.get(0);	
			int d2 = 0;
			for(int x3=1;x3<s.size();++x3){					
				if(min>s.get(x3)){
					min = s.get(x3);
					d2  = x3;
				}					
			}	
				
			
			m[i]=t.get(d2);
			s.clear();
			t.clear();
			
		}
	}

//This checks if previous k ie. tempk and current k are same.Used as terminating case.
	static int VerifyEqual(){	
		for(int i=0;i<p;++i){
			for(int j=0;j<n;++j){
				if(tempk[i][j]!=k[i][j]){
					return 0;
				}
			}				
		}		
		return 1;
	}	

	public static void main(String[] args) throws Exception {
		
		// TODO Auto-generated method stub
		Clustering();
		// C:\Users\buvan.suji\workspace\Distance\SampleData\results\TR_Fasta\OutputReads\chr11119005_1121390.fasta
		// C:\Users\buvan.suji\workspace\Distance\SampleData\results\TR_Fasta\OutputReads\chr1842862036_42864557.fasta
		//C:/Users/Subashchandran/Marseclipse/workspace/Data/TrResults/Finished Cluster/chr11119005_1121390.fasta
		//C:/Users/Subashchandran/Marseclipse/workspace/Data/TrResults/chr2440992_444078.fasta
		//C:/Users/Subashchandran/Marseclipse/workspace/Data/TrResults/chr68826819_8829365.fasta
	}



}
