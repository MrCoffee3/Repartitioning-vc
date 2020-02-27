package Util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class Group {
	public List<Vertex> newReplicaSets;
	public Set<Edge> edgesToSend;
	public int target;
	public int gain;

	public Group(List<Vertex> newReplicaSets, Set<Edge> edgesToSend) {
		this.newReplicaSets = newReplicaSets;
		this.edgesToSend = edgesToSend;
	}

	public Group(List<Vertex> newReplicaSets, Set<Edge> edgesToSend, int target, int gain) {
		this.newReplicaSets = newReplicaSets;
		this.edgesToSend = edgesToSend;
		this.target = target;
		this.gain = gain;
	}

	public static List<Vertex> getNewVerticesAfterSending(Set<Edge> tmpEdges, int exchangePartner, Partition p) {

		// construct new vertex replica sets of all vertices that change after edges
		// were sent
		HashMap<Integer, Vertex> newVertices = new HashMap<Integer, Vertex>();
		for (Edge e : tmpEdges) {
			List<Integer> edgeIDs = new ArrayList<Integer>();
			edgeIDs.add(e.u);
			edgeIDs.add(e.v);
			for (Integer id : edgeIDs) {
				if (!newVertices.containsKey(id)) {
					Vertex v = p.vertices.get(id);
					Vertex newV = new Vertex(v.vertexID, Util.copyList(v.replicaSet));
					boolean changed = false;
					if (p.edges.get(v.vertexID).size() == 0) {
						// after sending the edges, this vertex will be isolated and can be removed
						newV.replicaSet.remove(new Integer(p.id));
						changed = true;

					}
					if (!newV.replicaSet.contains(exchangePartner)) {
						newV.addReplica(exchangePartner);
						changed = true;
					}
					if (changed) {
						// System.out.println("Old replica set: " + v.replicaSet() + "\t" + v.master());
						// System.out.println("New replica set: " + newV.replicaSet() + "\t" +
						// newV.master());
						newVertices.put(newV.vertexID, newV);
					}
				}
			}
		}
		return new ArrayList<Vertex>(newVertices.values());
	}

	
}
