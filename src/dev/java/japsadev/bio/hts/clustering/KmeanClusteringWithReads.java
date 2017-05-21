package japsadev.bio.hts.clustering;

import java.util.ArrayList;
import java.io.*;

import japsa.seq.SequenceOutputStream;
import japsadev.bio.hts.clustering.PairDistance;
import japsadev.bio.hts.clustering.GettingTreadsFromFasta;

/**
 * @author buvan.suji
 *
 */

public class KmeanClusteringWithReads {	
	
	public static ArrayList<ArrayList<String>> 
		Clustering(ArrayList<String> reads) throws Exception{
		ArrayList<ArrayList<String>> tempClusterList = 
				new ArrayList<ArrayList<String>>();
		ArrayList<String> tempList = new ArrayList<String>();
			if(reads.size()<=2){
				tempClusterList.add(tempList);
				tempClusterList.add(tempList);
				tempClusterList.add(tempList);
				return tempClusterList;
			}
			else{
				int MaxReadLn = 0;
				int MinReadLn = 1000000;		
				
				int n;	
				int d[];
				int k[][];
				final int p =2;//number of clusters
				int tempk[][];
				int m[];
				int Nclusters = 0;
					
				double max = 0;
				double temp1 = 0;
				int index1=0;
				int index2=0;	
				int count1=0,count2=0;
				long startTime = System.nanoTime();		
				
				n = reads.size();					
				
				double[][] table = new double[n][n];	
				ArrayList<String> readsLengthRange = new ArrayList<String>();
				
				for (int x=0; x<n; ++x){
					for(int y=0;y<n;++y){
						table[x][y]=0;
					}
				}

				for (int x = 0; x < n; x++) {
					int temp = 0;
					for (int y = x + 1; y < n; y++) {				
						String r1 = reads.get(x);
						String r2 = reads.get(y);				
						PairDistance.EditDistanceResult result = PairDistance.compute(r1, r2);
						double normdist = ((double) result.getDistance())
								/ Math.max(r1.length(), r2.length());				
						table[x][y]=normdist; 
						temp = r1.length();
					}
					MaxReadLn = Math.max(MaxReadLn, temp);
					MinReadLn = Math.min(MinReadLn, temp);
				}
				readsLengthRange.add(""+MinReadLn);
				readsLengthRange.add(""+MaxReadLn);
							
				d=new int[n];
				
				 //name of the reads		
				for(int x=0;x<n;++x){
					d[x]=x;
				}	
				
				 //Initialising arrays 
				k=new int[p][n];
				tempk=new int[p][n];		
				m=new int [p];
					
				
				for (int x = 0; x < n; x++) {
					for (int y = x + 1; y < n; y++) {
						temp1 = table[x][y];
						if(max < temp1){
							max = temp1;
							index1=x;
							index2=y;
						}
					}
				}			
				
				 //Initializing m 
				m[0] = index1;
				m[1] = index2;		
				
				int temp=0;
				int flag=0;
				
				do{
					for(int x=0;x<p;++x){
						for(int y=0;y<n;++y){
							k[x][y]=-1;					
						}
					}
					
					for(int x=0;x<n;++x){
						temp=NewCluster(d[x], table, m);					
						if(temp==0){
							k[temp][count1++]=d[x];					
						}
						else if(temp==1){
							k[temp][count2++]=d[x];					
						}					 
					}
					
					NewMean(table, m, n, k); // call to method which will calculate mean at this step.
					
					
					flag = VerifyEqual(n, k, tempk); // check if terminating condition is satisfied.
					if(flag!=1){
						/*Take backup of k in tempk so that you can check for equivalence in next step*/
						for(int x=0;x<p;++x){
							for(int y=0;y<n;++y){
								tempk[x][y]=k[x][y];
							}
						}				
					}
					
					
					
					count1=0;count2=0;
					Nclusters += 1;
				}
				while(flag==0);
				long time = System.nanoTime()-startTime;
				double estimateTime = (double)time/1000000000.0;
				
				
				ArrayList<ArrayList<String>> clusterList = 
						new ArrayList<ArrayList<String>>();
				clusterList.add(readsLengthRange);
				
				for(int x=0;x<p;++x){
					ArrayList<String> tempcluster = new ArrayList<String>();				
					for(int y=0;k[x][y]!=-1 && y<n-1;++y){
						tempcluster.add(reads.get(k[x][y]));					
					}
					clusterList.add(tempcluster);				
				}
			return clusterList;	
			}	
			
	}
	

	
	public static int NewCluster(int a, double list[][], int m[] ){ 		
		
		double diff[] = new double[2];
			
		for(int i=0;i<2;++i){
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
		for(int i=0;i<2;++i){
			if(diff[i]<temp){
				temp=diff[i];
				val=i;
			}
			
		}//end of for loop
		
		return val;
	}
	
	static void NewMean(double list [][], int m[], int n, int k[][]){
		ArrayList<Integer> t = new ArrayList<Integer>() ;
		
		for(int i=0;i<2;++i){
			m[i]=0; // initializing means to 0			
		}
		
		//int cnt=0;
		for(int i=0;i<2;++i){			
			
			for(int j=0;j<n-1;++j){								
				if(k[i][j]!=-1){
					t.add(k[i][j]);
				}
			}
			if(t.size()>1){
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
			}else{
				m[i]=t.get(0);
			}	
			
		}
	}

//This checks if previous k ie. tempk and current k are same.Used as terminating case.
	static int VerifyEqual(int n, int k[][], int tempk[][]){	
		for(int i=0;i<2;++i){
			for(int j=0;j<n;++j){
				if(tempk[i][j]!=k[i][j]){
					return 0;
				}
			}				
		}		
		return 1;
	}	

	



}
