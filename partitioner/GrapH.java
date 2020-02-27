package partitioner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import Util.*;

public class GrapH {
	private Partition[] partitionList;
	// public Set<Vertex> vertexinC;
	// public Set<Edge> edgelist;
	private int pNum;
	private float balance;
	private int vertexNum;
	private int edgeNum;
	private float maxSize;
	private String edgePath;
	private int migration = 0;
	private int currentMigration = 0;
	private int currentSearched = 0;
	private int currentfound = 0;
	private List<HashMap<Integer, Integer>> backoffUntilIterations;

	// private HashMap<Integer, Integer> backoffUntilIteration;
	private int count = 0;
	private double RF;
	private CyclicBarrier cb1;
	private CyclicBarrier cb2;

	public GrapH(int pNum, float balance, int vNum, String iEP) {
		// TODO Auto-generated constructor stub
		this.pNum = pNum;
		this.balance = balance;
		this.vertexNum = vNum;
		this.edgePath = iEP;
		partitionList = new Partition[pNum];
		this.backoffUntilIterations = new ArrayList<HashMap<Integer, Integer>>();
		for (int i = 0; i < pNum; i++) {
			partitionList[i] = new Partition(i);
			backoffUntilIterations.add(new HashMap<Integer, Integer>());
			for (int j = 0; j < pNum; j++) {
				if (i != j) {
					this.backoffUntilIterations.get(i).put(j, 0);
				}
			}
		}
		// this.backoffUntilIteration = new HashMap<Integer, Integer>();
		// for (int i = 0; i < pNum; i++) {
		// this.backoffUntilIteration.put(i, 0);
		// }
	}

	public void load() throws IOException {
		// String str1;
		Set<Integer> set = new HashSet<Integer>();
		BufferedReader in1;// 边
		for (int i = 0; i < pNum; i++) {
			String fileName = "partition" + i + ".txt";
			in1 = new BufferedReader(new FileReader(new File(edgePath + "\\" + fileName)));// 初始边
			String string;
			while ((string = in1.readLine()) != null) {
				edgeNum++;
				String[] vertexs = string.split(" ");
				int u = Integer.parseInt(vertexs[0]);
				int v = Integer.parseInt(vertexs[1]);
				set.add(u);
				set.add(v);
				if (!partitionList[i].vertices.containsKey(u)) {
					partitionList[i].newVertex(u);
				}
				if (!partitionList[i].vertices.containsKey(v)) {
					partitionList[i].newVertex(v);
				}
				Edge edge = new Edge(u, v);
				partitionList[i].addEdgeToPartition(edge);
			}
			in1.close();
		}
		vertexNum = set.size();
		maxSize = (float) (edgeNum * 1.0 / pNum * balance);
		exchangeRepSet();
	}

	public void exchangeRepSet() {
		int totalVnum = 0;
		for (Partition send : partitionList) {
			totalVnum += send.vertices.size();
			for (Partition re : partitionList) {
				if (re.id != send.id) {
					for (Map.Entry<Integer, Vertex> entry : send.vertices.entrySet()) {
						int key = entry.getKey();
						if (re.vertices.containsKey(key)) {
							Vertex v = re.vertices.get(key);
							if (!v.replicaSet.contains(send.id)) {
								v.addReplica(send.id);
							}
						}
					}
				}
			}
		}
		RF = totalVnum * 1.0 / vertexNum;
		System.out.println("RF:" + totalVnum * 1.0 / vertexNum);
	}

	public int getPartner(Partition p) {
		List<Integer> candidates = new ArrayList<Integer>();
		for (int i = 0; i < pNum; i++) {
			if (i != p.id) {
				candidates.add(i);
			}
		}
		/*
		 * Remove candidates that have no exchanged traffic, because no improvement can
		 * be found. Remove also candidates that are currently ignored.
		 */
		PartnerComparator comparator = new PartnerComparator(p, pNum);

		/*
		 * Remove candidates that have no exchanged traffic, because no improvement can
		 * be found. Remove also candidates that are currently ignored.
		 */
		Iterator<Integer> iter = candidates.iterator();
		while (iter.hasNext()) {
			Integer tmpPartner = iter.next();
			if (comparator.exchangedTraffic(tmpPartner) <= 0.1) {
				iter.remove();
			}
			if (count < backoffUntilIterations.get(p.id).get(tmpPartner)) {
				iter.remove();
			}
			/*
			 * if(partitionList[tmpPartner].w<p.w) { iter.remove(); }
			 */
		}

		/*
		 * Sort candidates by exchanged traffic
		 */
		Collections.sort(candidates, comparator);
		// System.out.println("Partner candidates (sorted): " + candidates);
		if (candidates.size() == 0) {
			return -1;
		} else {
			backoffUntilIterations.get(p.id).put(candidates.get(0), count + 1);
			// System.out.println(p.id+"->"+candidates.get(0));
			return candidates.get(0);
		}
	}

	public boolean hasLock(Vertex v, Partition p) {
		// if (p.isMaster(v)) {
		/*
		 * Default: one partition has all locks it needs
		 */
		int partitionHasLock = count % pNum;
		if (v.replicaSet.contains(partitionHasLock)) {
			if (p.id == partitionHasLock) {
				return true;
			} else {
				return false;
			}
		} else {
			/*
			 * Default partition does not need this lock. Give it to the partition in the
			 * replica set with
			 */
			List<Integer> replicaSet = new ArrayList<Integer>(v.replicaSet);
			partitionHasLock = replicaSet.get(count % replicaSet.size());
			if (partitionHasLock == p.id) {
				return true;
			} else {
				return false;
			}
		}
	}

	private Group initExchangeCandidates(Partition p) {

		int exchangePartner = getPartner(p);
		// Consider only adjacent vertices
		List<Vertex> vertexCandidates = new ArrayList<Vertex>();
		// for (Vertex v : p.subgraph().getVertices().values()) {
		for (Iterator<Map.Entry<Integer, Vertex>> it = p.vertices.entrySet().iterator(); it.hasNext();) {
			Vertex v = it.next().getValue();
			if (v.replicaSet.contains(exchangePartner) && hasLock(v, p)) {
				vertexCandidates.add(v);
			}
		}

		List<Group> bagList = new ArrayList<Group>();
		if (vertexCandidates.size() > 0) {
			Collections.sort(vertexCandidates, new Comparator<Vertex>() {
				@Override
				public int compare(Vertex o1, Vertex o2) {
					// TODO Auto-generated method stub
					if (o1.replicaSet.size() > o2.replicaSet.size()) {
						return -1;
					} else if (o1.replicaSet.size() < o2.replicaSet.size()) {
						return 1;
					} else {
						return 0;
					}
				}
			});
			double cap = maxSize - partitionList[exchangePartner].w;
			for (Vertex v : vertexCandidates) {
				if (cap < 0) {
					break;
				}
				Set<Edge> tmpEdges = new HashSet<Edge>();

				/*
				 * Add all in edges
				 */
				Iterator<Integer> iter2 = p.edges.get(v.vertexID).iterator();
				while (iter2.hasNext()) {
					Edge edge = new Edge(iter2.next(), v.vertexID);
					tmpEdges.add(edge);
					currentSearched++;
				}

				// consider only edges where lock on both endpoints are held
				Iterator<Edge> iter = tmpEdges.iterator();
				while (iter.hasNext()) {
					Edge e = iter.next();
					if (!hasLock(p.vertices.get(e.u), p) || !hasLock(p.vertices.get(e.v), p)) {
						iter.remove();
					}
				}

				// remove all edges to send
				for (Edge e : tmpEdges) {
					p.removeEdge(e);
				}
				List<Vertex> tmpVertices = Group.getNewVerticesAfterSending(tmpEdges, exchangePartner, p);
				int evc = 0;
				int ivc = 0;
				for (Vertex newV : tmpVertices) {
					Vertex oldV = p.vertices.get(newV.vertexID);
					if (oldV.replicaSet.contains(exchangePartner)) {
						evc++;
					}
					if (newV.replicaSet.contains(p.id)) {
						ivc++;
					}
				}

				if (evc - ivc > 0) {
					// System.out.println(evc-ivc);
					cap -= tmpEdges.size();
					Group bag = new Group(tmpVertices, tmpEdges);
					bagList.add(bag);
				} else {
					for (Edge e : tmpEdges) {
						p.addEdgeToPartition(e);
					}
				}
			}
		}

		// Now we have a list of bag of edges, which individual exchange results in
		// lower costs
		// if we send the union, this results in lower costs, too.
		Set<Edge> unionEdges = new HashSet<Edge>();
		List<Vertex> unionVertices;
		for (Group bag : bagList) {
			for (Edge e : bag.edgesToSend) {
				unionEdges.add(e);
			}
			// for (Vertex v : bag.newReplicaSets) {
			// // the newer vertex version replaces the older one.
			// unionVertices.put(v.id(), v);
			// }
		}
		unionVertices = Group.getNewVerticesAfterSending(unionEdges, exchangePartner, p);
		int evc = 0;
		int ivc = 0;
		for (Vertex newV : unionVertices) {
			Vertex oldV = p.vertices.get(newV.vertexID);
			if (oldV.replicaSet.contains(exchangePartner)) {
				evc++;
			}
			if (newV.replicaSet.contains(p.id)) {
				ivc++;
			}
		}
		Group uniongroup = new Group(unionVertices, unionEdges, exchangePartner, evc - ivc);
		return uniongroup;

	}

	public void repartitioning() {
		long t1 = System.currentTimeMillis();
		double currentRF = RF;
		ExecutorService pool = Executors.newCachedThreadPool();
		cb1 = new CyclicBarrier(pNum);
		cb2 = new CyclicBarrier(pNum + 1);
		while (count == 0 || currentRF < RF) {
			count++;
			// List<repartitionThread> list=new ArrayList<repartitionThread>();
			for (Partition partition : partitionList) {
				// initExchangeCandidates(partition);
				// repartitionThread rThread=new repartitionThread(partition);
				// list.add(rThread);
				// rThread.start();
				pool.execute(new repartitionThread(partition));
			}
			try {
				cb2.await();
			} catch (InterruptedException | BrokenBarrierException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			/*
			 * for (repartitionThread rThread : list) { try { rThread.join(); } catch
			 * (InterruptedException e) { // TODO Auto-generated catch block
			 * e.printStackTrace(); } }
			 */
			// System.out.println(count+": "+migrationCurrent);
			// System.out.println(count+"finish");
			int totalVnum = 0;
			for (Partition send : partitionList) {
				totalVnum += send.vertices.size();
			}
			RF = currentRF;
			currentRF = totalVnum * 1.0 / vertexNum;
			// System.out.println("superstep:" + count + "\tRF:" + currentRF + "\t" +
			// currentSearched + "\t" + currentfound
			// + "\t" + currentMigration);
			currentMigration = 0;
			currentfound = 0;
			currentSearched = 0;

		}
		long t2 = System.currentTimeMillis();
		int totalVnum = 0;
		for (Partition send : partitionList) {
			totalVnum += send.vertices.size();
		}
		RF = currentRF;
		currentRF = totalVnum * 1.0 / vertexNum;
		System.out.println("iterations:" + count);
		System.out.println("RF:" + totalVnum * 1.0 / vertexNum);
		float maxw = 0;
		for (Partition partition : partitionList) {
			if (partition.w > maxw) {
				maxw = partition.w;
			}
		}
		float avg = (float) edgeNum / (float) pNum;
		System.out.println("maxbal:" + ((maxw - avg) / avg));
		System.out.println("migration:" + migration);
		System.out.println("run time:" + (t2 - t1) * 1.0 / 1000 + "s");
	}

	class repartitionThread extends Thread {
		private Partition p;

		public repartitionThread(Partition p) {
			// TODO Auto-generated constructor stub
			this.p = p;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			super.run();
			Group group = initExchangeCandidates(p);
			currentfound += group.edgesToSend.size();
			try {
				// System.out.println("partition "+p.id+"finish finding group");
				cb1.await();
			} catch (InterruptedException | BrokenBarrierException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if (group.gain > 0) {
				for (Vertex newV : group.newReplicaSets) {
					Vertex oldV = p.vertices.get(newV.vertexID);
					if (!oldV.replicaSet.contains(group.target)) {
						partitionList[group.target].vertices.put(newV.vertexID, newV);
					}
					for (int pid : newV.replicaSet) {
						Vertex vertex = partitionList[pid].vertices.get(newV.vertexID);
						if (vertex != null) {
							vertex.replicaSet = Util.copyList(newV.replicaSet);
						}
					}

					if (!newV.replicaSet.contains(p.id)) {
						p.removeVertex(newV.vertexID);
					}
				}
				for (Edge edge : group.edgesToSend) {
					partitionList[group.target].addEdgeToPartition(edge);
					migration++;
					currentMigration++;
				}
			}
			try {
				// System.out.println("partition "+p.id+"finish migrating group");
				cb2.await();
			} catch (InterruptedException | BrokenBarrierException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String[] datas = { /*"dblp", "roadNet-PA", "youtube", "lj",*/ "orkut" };
		int[] ps = { /* 2, 4, 8, */16 };
		String[] trackers = { "RBSEP"/*, "random"*/ };
		String[] initialMehods = { /*"Ne",*/ "DBH" };
		for (String data : datas) {
			for (String initial : initialMehods) {
				for (String streaming : trackers) {
					for (int i : ps) {
						int pNum = i;
						float balance = 1.1f;// 平衡系数
						int vNum = 0;// 点数量317080,1088092,1134890,3997962,7600000,3072441
						switch (data) {
						case "dblp":
							vNum = 317080;
							break;
						case "roadNet-PA":
							vNum = 1088092;
							break;
						case "youtube":
							vNum = 1134890;
							break;
						case "lj":
							vNum = 3997962;
							break;
						case "sp":
							vNum = 7600000;
							break;
						case "orkut":
							vNum = 3072441;
							break;
						default:
							break;
						}
						String iEP = "";

						iEP = "D:\\dataset\\Initial+dynamic\\" + initial + "+" + streaming + "\\" + data
								+ "2-partition-" + pNum;
						System.out.println("-----------------" + data + "-" + pNum + "-" + initial + "+" + streaming
								+ "-----------------");
						GrapH partitioner = new GrapH(pNum, balance, vNum, iEP);
						try {
							partitioner.load();
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						partitioner.repartitioning();
					}
				}
			}
		}
	}

}
