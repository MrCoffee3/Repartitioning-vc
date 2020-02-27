package Util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class Partition {
	public int id;
	public HashMap<Integer, Vertex> vertices = new HashMap<Integer, Vertex>(100000); // TODO: change to trove
	public HashMap<Integer, HashSet<Integer>> edges = new HashMap<Integer, HashSet<Integer>>(100000);
	// public HashMap<Integer, Double> dis;
	public int w = 0;

	public Partition(int i) {
		// TODO Auto-generated constructor stub
		this.id = i;
	}

	public void addEdgeToPartition(Edge e) {
		if (!edges.containsKey(e.v)) {
			edges.put(e.v, new HashSet<Integer>());
		}
		if (!edges.containsKey(e.u)) {
			edges.put(e.u, new HashSet<Integer>());
		}
		if (!edges.get(e.v).contains(e.u)) {
			edges.get(e.v).add(e.u);
			edges.get(e.u).add(e.v);
			w++;
		}
	}

	public void newVertex(int v_id) {
		List<Integer> l = new ArrayList<Integer>();
		l.add(id);
		Vertex v = new Vertex(v_id, l);
		vertices.put(v_id, v);
		edges.put(v_id, new HashSet<Integer>());
		edges.put(v_id, new HashSet<Integer>());
	}

	public void removeEdge(Edge e) { // TODO: input id1, id2 instead of wrapper edge
		if (edges.containsKey(e.v) && edges.containsKey(e.u)) {
			edges.get(e.v).remove(e.u);
			edges.get(e.u).remove(e.v);
			w--;
		}
	}

	public void removeVertex(Integer v_id) {
		// TODO: remove debug

		if (!vertices.containsKey(v_id) || !edges.containsKey(v_id)) {
			System.out.println("Vertex not present: Error!! " + v_id);
		}
		vertices.remove(v_id);
		edges.remove(v_id);

	}
}
