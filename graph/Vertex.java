package Util;

import java.util.List;
import java.util.Set;

public class Vertex {
	public int vertexID;
	public int dis;
	public int notInGroup;
	public List<Integer> replicaSet;
	public int checked=0;
	public Vertex(int id, List<Integer> replicaSet) {
		// TODO Auto-generated constructor stub
		this.vertexID=id;
		this.replicaSet=replicaSet;
		this.dis=10000;
	}
	
	public void addReplica(Integer replica) {
		int pos = 0;
		for (pos=0; pos<replicaSet.size(); pos++) {
			if (replicaSet.get(pos)>replica) {
				break;
			}
		}
		replicaSet.add(pos, replica);
	}
	
	public boolean allInGroup(Set<Integer> set, Partition p) {
		for(int i : p.edges.get(vertexID)) {
			if(!set.contains(i)) {
				return false;
			}
		}
		return true;
	}
}
