package Util;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;


public class PartnerComparator implements Comparator<Integer> {

	private HashMap<Integer, Integer> commonReplicas;
	private HashMap<Integer, Double> exchangedTraffic;
	
	public PartnerComparator(Partition p, int pNum) {
		this.commonReplicas = new HashMap<Integer, Integer>();
		this.exchangedTraffic = new HashMap<Integer, Double>();
		for (int partitionID=0;partitionID<pNum;partitionID++) {
			this.commonReplicas.put(partitionID, 0);
			this.exchangedTraffic.put(partitionID, 0.0d);
		}
//		for (Vertex v : p.subgraph().getVertices().values()) {
		for ( Iterator<Map.Entry<Integer, Vertex>> it = p.vertices.entrySet().iterator(); it.hasNext(); ) {
			Vertex v = it.next().getValue();
			for (Integer re : v.replicaSet) {
				int newVal = this.commonReplicas.get(re) + 1;
				this.commonReplicas.put(re, newVal);
				
				double newTraffic = this.exchangedTraffic.get(re) + 1;
				this.exchangedTraffic.put(re, newTraffic);
			}
		}
		//System.out.println("Common replicas: " + commonReplicas);
		//System.out.println("Exchanged traffic: " + exchangedTraffic);
		/*for (Map.Entry<Integer, Double> entry : exchangedTraffic.entrySet()) {
			System.out.println(entry.getKey()+";"+entry.getValue());
		}*/
	}
	
	/**
	 * Returns the number of common replicas this partition
	 * has with partition p
	 * @param p
	 * @return
	 */
	public int commonReplicas(Integer p) {
		return commonReplicas.get(p);
	}
	
	/**
	 * Returns the exchanged traffic between this partition and p
	 * @param p
	 * @return
	 */
	public double exchangedTraffic(Integer p) {
		return exchangedTraffic.get(p);
	}

	
	@Override
	public int compare(Integer arg0, Integer arg1) {
		double traffic0 = exchangedTraffic.get(arg0);
		double traffic1 = exchangedTraffic.get(arg1);
		
		if (traffic0<traffic1) {
			return 1;
		} else if (traffic0==traffic1) {
			return 0;
		} else {
			return -1;
		}
	}

	
}
