package japsa.bio.tr;

import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.stat.clustering.Cluster;
import org.apache.commons.math3.stat.clustering.Clusterable;


public class ReadAllele implements Clusterable<ReadAllele> {
	String readName = "";
	public double copy_number;


	public ReadAllele(String name, double cn){
		readName = name + "";//copy
		copy_number = cn;
	}

	public ReadAllele(double cn){
		copy_number = cn;
	}

	public double distanceFrom(ReadAllele p) {
		return Math.abs(copy_number - p.copy_number);
	}

	public ReadAllele centroidOf(Collection<ReadAllele> pp) {
		double sum = 0.0;
		for (ReadAllele p:pp){
			sum += p.copy_number;
		}
		return new ReadAllele(sum/pp.size());
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof ReadAllele)) {
			return false;
		}
		ReadAllele p = (ReadAllele) o;

		return copy_number == p.copy_number;
	}

	public static double inertia(List<Cluster<ReadAllele>> clusters){
		double sum = 0.0;
		int count = 0;
		for (Cluster<ReadAllele> cluster: clusters){
			ReadAllele cenroid = cluster.getCenter();
			for (ReadAllele allele:cluster.getPoints()){
				sum += allele.distanceFrom(cenroid);
				count ++;
			}
		}
		return sum / count;
	}
}