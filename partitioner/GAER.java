package partitioner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import Util.Edge;
import Util.Partition;
import Util.Util;
import Util.Vertex;
import Util.Group;

public class GAER {
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
	// private HashMap<Integer, Integer> backoffUntilIteration;
	private int count = 0;
	private int maxGroupSize;
	public float T0 = 0.4f;
	public float d = 0.03f;
	public float Tr = 1;
	private double RF;
	private int checked;
	private float gainRF;
	int notMigreteByW = 0;
	int currentSearched = 0;

	public HashMap<Integer, Queue<Group>> AllgroupList = new HashMap<Integer, Queue<Group>>();
	private CyclicBarrier cb;
	private int currentfound = 0;

	public GAER(int pNum, float balance, int vNum, String iEP, int mgz, int checked, float gainRF) {
		// TODO Auto-generated constructor stub
		this.pNum = pNum;
		this.balance = balance;
		this.vertexNum = vNum;
		this.edgePath = iEP;
		this.maxGroupSize = mgz;
		partitionList = new Partition[pNum];
		this.checked = checked;
		this.gainRF = gainRF;
		for (int i = 0; i < pNum; i++) {
			partitionList[i] = new Partition(i);
		}

		// this.backoffUntilIteration = new HashMap<Integer, Integer>();
		// for (int i = 0; i < pNum; i++) {
		// this.backoffUntilIteration.put(i, 0);
		// }
	}

	public void load() throws IOException {
		// String str1;
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

	public List<Vertex> markDis(Partition partition) {
		// partition.dis = new HashMap<Integer, Double>();
		List<Vertex> vertexs = new ArrayList<Vertex>();
		Queue<Integer> queue = new LinkedList<Integer>();
		// Set<Vertex> vertexs=new HashSet<Vertex>();
		for (Integer vertex : partition.vertices.keySet()) {
			if (partition.vertices.get(vertex).replicaSet.size() > 1) {
				vertexs.add(partition.vertices.get(vertex));
				queue.add(vertex);
				// partition.dis.put(vertex, -1d);
			}
		}
		int dis = 0;
		while (!queue.isEmpty()) {
			int s = queue.size();
			// System.out.println(s);
			for (int i = 0; i < s; i++) {
				int vertex = queue.poll();
				// System.out.println("error dis!"+partition.id+";"+dis);
				partition.vertices.get(vertex).dis = dis;
				partition.vertices.get(vertex).notInGroup = partition.edges.get(vertex).size();
				// partition.dis.put(vertex, (1.0 / (dis + 1)));
				for (Iterator<Integer> it = partition.edges.get(vertex).iterator(); it.hasNext();) {
					int neighbor = it.next();
					if (partition.vertices.get(neighbor).dis > dis) {
						queue.add(neighbor);
						// partition.dis.put(neighbor, -1d);
					}
				}
			}
			// System.out.println(queue.size()+";"+partition.dis.size());
			dis++;
		}
		/*
		 * if (partition.dis.size() != partition.vertices.size()) {
		 * System.out.println(partition.id + ";" + count + ";" + partition.dis.size() +
		 * ";" + partition.vertices.size() + ";markDis error!!"); System.exit(-1); }
		 */
		return vertexs;
	}

	public Group findGroup(Vertex vertex, Partition p) {
		int[] caps = new int[pNum];
		for (int i = 0; i < pNum; i++) {
			caps[i] = (int) maxSize - partitionList[i].w;
		}
		HashMap<Integer, Double> scoreMap = new HashMap<Integer, Double>();
		// Set<Integer> needReset = new HashSet<Integer>();
		Set<Integer> vList = new HashSet<Integer>();
		Set<Edge> edgeGroup = new HashSet<Edge>();
		PriorityQueue<Integer> queue = new PriorityQueue<Integer>(new Comparator<Integer>() {

			@Override
			public int compare(Integer o1, Integer o2) {
				// TODO Auto-generated method stub
				if (scoreMap.get(o1) > scoreMap.get(o2)) {
					return -1;
				} else if (scoreMap.get(o1) < scoreMap.get(o2)) {
					return -1;
				}
				return 0;
			}
		});

		int[] ovc = new int[pNum];
		int maxovc = 0;
		int maxp = -1;
		int ivc = 0;
		int groupsize = 0;
		scoreMap.put(vertex.vertexID, 1.0 / (vertex.dis + 1));
		queue.add(vertex.vertexID);
		currentSearched++;

		while (!queue.isEmpty()) {
			int vid = queue.poll();
			/*
			 * if(vList.contains(vid)) { System.out.println("?"); continue; }
			 */
			Vertex v = p.vertices.get(vid);
			HashSet<Vertex> edges = new HashSet<Vertex>();
			for (Integer n : p.edges.get(vid)) {
				Vertex nVertex = p.vertices.get(n);
				// if(nVertex.dis!=-1&&(nVertex.dis!=10000||v.dis!=10000)) {
				edges.add(nVertex);
				// }
			}
			vList.add(vid);
			if (vList.size() == 1) {
				for (int i : v.replicaSet) {
					if (i != p.id) {
						ovc[i]++;
					}
				}
				if (edges.size() > 0) {
					ivc++;
				}
				/*
				 * int highest=2; int[] highestVertex = new int[highest]; double[] highScore =
				 * new double[highest];
				 */
				for (Iterator<Vertex> it = edges.iterator(); it.hasNext();) {
					Vertex neigh = it.next();
					/*
					 * if (neigh.replicaSet.size() > 0 && !hasLock(neigh, p)) { continue; }
					 */
					if (neigh.dis != 10000) {
						// needReset.add(neighbor);
						neigh.notInGroup--;
						double newScore = 1.0 / (neigh.dis + 1) + 1.0 / v.notInGroup
								+ 1.0 / Math.pow(neigh.notInGroup + 1, 2);
						// if (newScore > highScore[0]) {
						// highestVertex[1] = highestVertex[0];
						// highScore[1] = highScore[0];
						// highestVertex[0] = neigh.vertexID;
						// highScore[0] = newScore;
						// } else if (newScore > highScore[1]) {
						// highestVertex[1] = neigh.vertexID;
						// highScore[1] = newScore;
						// }
						scoreMap.put(neigh.vertexID, newScore);

						queue.add(neigh.vertexID);
						currentSearched++;
					}

				}
				/*
				 * for(int k=0;k<highest;k++) { if(highestVertex[k]!=0) {
				 * scoreMap.put(highestVertex[k], highScore[k]); queue.offer(highestVertex[k]);
				 * currentSearched ++; } }
				 */

			} else {

				int eNum = 0;
				int ivcReduce = 0;
				// Set<Edge> tmpEdges = new HashSet<Edge>();
				for (Vertex neigh : edges) {
					/*
					 * if (neigh.replicaSet.size() > 0 && !hasLock(neigh, p)) { continue; }
					 */
					if (neigh.dis != 10000) {
						// needReset.add(i);
						neigh.notInGroup--;
						scoreMap.put(neigh.vertexID, 0d);
						if (vList.contains(neigh.vertexID)) {
							Edge edge = new Edge(neigh.vertexID, vid);
							edgeGroup.add(edge);
							// tmpEdges.add(edge);
							// p.removeEdge(edge);
							eNum++;
							if (neigh.notInGroup == 0) {

								ivcReduce++;
							}
						}
					}
				}
				/*
				 * for (Edge edge : tmpEdges) { p.removeEdge(edge); }
				 */

				groupsize += eNum;
				if (groupsize <= maxGroupSize) {
					if (v.replicaSet.size() > 1) {
						for (int i : v.replicaSet) {
							if (i != p.id) {
								ovc[i]++;
								if (ovc[i] > maxovc && caps[i] > 0) {
									maxovc = ovc[i];
									maxp = i;
								}
							}
						}
					}
					if (eNum < edges.size()) {
						ivc++;
					}
					ivc -= ivcReduce;
					if (maxovc * Tr > ivc && maxp != -1) {
						List<Vertex> tmpVertices = new ArrayList<Vertex>();
						for (Integer vInteger : vList) {
							tmpVertices.add(p.vertices.get(vInteger));
						}
						/*
						 * List<Vertex> t = Group.getNewVerticesAfterSending(edgeGroup, maxp, p);
						 * if(t.size()!=tmpVertices.size()) { for (Vertex vertex2 : t) {
						 * System.out.print(vertex2.vertexID+" "); } System.out.println(); for (Integer
						 * vertex2 : vList) { System.out.print(vertex2+" "); } System.out.println(); for
						 * (Edge edge : edgeGroup) { System.out.print(edge.v+" "+edge.u+";"); }
						 * System.out.println(); }
						 */
						for (Vertex v1 : tmpVertices) {
							v1.dis = 10000;
							/*
							 * for (int i : p.edges.get(v1.vertexID)) { if(!vList.contains(i)) {
							 * v1.dis=10000; break; } }
							 */
						}
						for (int v1 : scoreMap.keySet()) {
							if (!vList.contains(v1)) {
								reset(v1, p);
							}
						}
						caps[maxp] -= edgeGroup.size();
						currentfound += edgeGroup.size();

						int evc = 0;
						int vc = 0;
						for (Vertex vertex2 : tmpVertices) {
							if (vertex2.replicaSet.contains(maxp)) {
								evc++;
							}
							for (int i : p.edges.get(vertex2.vertexID)) {
								if (!tmpVertices.contains(p.vertices.get(i))) {
									vc++;
									break;
								}
							}
						}
						if (evc != maxovc) {
							System.out.println("evc");
						}
						if (vc != ivc) {
							System.out.println(ivc + ";" + vc);
							for (Vertex vertex2 : tmpVertices) {
								System.out.print(vertex2.vertexID + " " + p.edges.get(vertex2.vertexID).size() + " "
										+ p.id + ";");
							}
							System.out.println();
							for (Edge edge : edgeGroup) {
								System.out.print(edge.v + " " + edge.u + ";");
							}
							System.out.println();
						}

						/*
						 * for (Edge edge : edgeGroup) { p.removeEdge(edge); }
						 */
						return new Group(tmpVertices, edgeGroup, maxp, maxovc - ivc);
					} else {
						/*
						 * int highest=2; int[] highestVertex = new int[highest]; double[] highScore =
						 * new double[highest];
						 */
						for (Iterator<Vertex> it = edges.iterator(); it.hasNext();) {
							Vertex neighbor = it.next();
							/*
							 * if (neighbor.replicaSet.size() > 0 && !hasLock(neighbor, p)) { continue; }
							 */
							if (neighbor.dis != 10000 && !vList.contains(neighbor.vertexID)) {
								double score1 = 1.0 / (neighbor.dis + 1);
								double score2 = 0;
								for (int neiborInGroup : p.edges.get(neighbor.vertexID)) {
									if (vList.contains(neiborInGroup)) {

										score2 += 1.0 / p.vertices.get(neiborInGroup).notInGroup;
									}
								}
								double score3 = 1.0 / Math.pow(neighbor.notInGroup + 1, 2);
								double newScore = score1 + score2 + score3;

								/*
								 * if (newScore > highScore[0]) { highestVertex[1] = highestVertex[0];
								 * highScore[1] = highScore[0]; highestVertex[0] = neighbor.vertexID;
								 * highScore[0] = newScore; } else if (newScore > highScore[1]) {
								 * highestVertex[1] = neighbor.vertexID; highScore[1] = newScore; }
								 */

								scoreMap.put(neighbor.vertexID, newScore);
								queue.remove(neighbor.vertexID);
								queue.add(neighbor.vertexID);
								currentSearched++;
							}
						}
						/*
						 * for(int k=0;k<highest;k++) { if(highestVertex[k]!=0) {
						 * scoreMap.put(highestVertex[k], highScore[k]); queue.remove(highestVertex[k]);
						 * queue.offer(highestVertex[k]); currentSearched ++; } }
						 */
					}
				} else {
					/*
					 * for (Edge edge : edgeGroup) { p.addEdgeToPartition(edge); }
					 */
					for (int v1 : scoreMap.keySet()) {
						reset(v1, p);
					}
					return null;
				}
			}
		}
		/*
		 * for (Edge edge : edgeGroup) { p.addEdgeToPartition(edge); }
		 */
		for (int v : scoreMap.keySet()) {
			reset(v, p);
		}
		return null;

	}

	/*
	 * public boolean hasLock(Vertex v, Partition p) { // if (p.isMaster(v)) {
	 * 
	 * Default: one partition has all locks it needs
	 * 
	 * int partitionHasLock = count % pNum; if
	 * (v.replicaSet.contains(partitionHasLock)) { if (p.id == partitionHasLock) {
	 * return true; } else { return false; } } else {
	 * 
	 * Default partition does not need this lock. Give it to the partition in the
	 * replica set with
	 * 
	 * List<Integer> replicaSet = new ArrayList<Integer>(v.replicaSet);
	 * partitionHasLock = replicaSet.get(count % replicaSet.size()); if
	 * (partitionHasLock == p.id) { return true; } else { return false; } } }
	 */

	public void repartitioning() throws InterruptedException, ExecutionException {
		/*
		 * List<repartitionThread> list=new ArrayList<repartitionThread>(); for
		 * (Partition p : partitionList) { repartitionThread rThread=new
		 * repartitionThread(p); rThread.setName("Partiton"+p.id); list.add(rThread); }
		 */
		ExecutorService pool = Executors.newCachedThreadPool();
		cb = new CyclicBarrier(pNum + 1);
		long t1 = System.currentTimeMillis();

		double currentRF = RF;

		while (count == 0 || Tr < 1 || RF - currentRF > gainRF) {
			cb.reset();
			Tr = Math.min(1, (float) (Math.round((T0 + count * d) * 100)) / 100);
			count++;
			long ts = System.currentTimeMillis();
			for (int i1 = 0; i1 < pNum; i1++) {
				repartitionThread rThread = new repartitionThread(partitionList[i1]);
				pool.submit(rThread);
			}
			try {
				cb.await();
			} catch (BrokenBarrierException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			long te = System.currentTimeMillis();
			int currentMigration = 0;
			for (Partition p : partitionList) {
				Queue<Group> groups = AllgroupList.get(p.id);
				while (!groups.isEmpty()) {
					Group group = groups.poll();
					int k = group.target;
					if (group.target == p.id) {
						System.out.println("error in findgroup: migrate to itself  ");
					}

					for (Edge edge : group.edgesToSend) {
						p.removeEdge(edge);
					}
					List<Vertex> tmpVertices = Group.getNewVerticesAfterSending(group.edgesToSend, group.target, p);
					group.newReplicaSets = tmpVertices;
					if (checkGroup(group, p)) {
						currentMigration += group.edgesToSend.size();
						for (Vertex newV : group.newReplicaSets) {
							Vertex oldV = p.vertices.get(newV.vertexID);
							if (!oldV.replicaSet.contains(k)) {
								partitionList[k].vertices.put(newV.vertexID, newV);
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
							partitionList[k].addEdgeToPartition(edge);
							migration++;
						}
					} else {
						notMigreteByW += group.edgesToSend.size();
						for (Edge edge : group.edgesToSend) {
							p.addEdgeToPartition(edge);
						}
					}

				}

			}
			AllgroupList.clear();
			int totalVnum = 0;
			for (Partition send : partitionList) {
				totalVnum += send.vertices.size();
			}
			RF = currentRF;
			currentRF = totalVnum * 1.0 / vertexNum;
			System.out.println("superstep:" + count + "\tRF:" + totalVnum * 1.0 / vertexNum + "\t" + currentSearched
					+ "\t" + currentfound + "\t" + currentMigration + "\t" + notMigreteByW + "\tTr:" + Tr + "\tmgz:"
					+ maxGroupSize + "\t" + (te - ts) * 1.0 / 1000);
			currentfound = 0;
			currentSearched = 0;
			notMigreteByW = 0;
			if (maxGroupSize < 200)
				maxGroupSize += 5;

		}
		long t2 = System.currentTimeMillis();
		int totalVnum = 0;
		for (Partition send : partitionList) {
			totalVnum += send.vertices.size();
		}
		System.out.println("Total iterations:" + count);
		System.out.println("RF:" + totalVnum * 1.0 / vertexNum);
		System.out.println("migration:" + migration);
		int totalw = 0;
		for (Partition partition : partitionList) {
			totalw += partition.w;
		}
		System.out.println("edgeNum:" + totalw + "\tvertexNum:" + vertexNum);
		System.out.println("run time:" + (t2 - t1) * 1.0 / 1000 + "s");
	}

	public boolean checkGroup(Group group, Partition belong) {
		if (partitionList[group.target].w > maxSize) {
			return false;
		}
		int evc = 0;
		int ivc = 0;
		for (Vertex newV : group.newReplicaSets) {
			Vertex oldV = belong.vertices.get(newV.vertexID);
			if (oldV.replicaSet.contains(group.target)) {
				evc++;
			}
			if (newV.replicaSet.contains(belong.id)) {
				ivc++;
			}
		}
		return evc > ivc;
	}

	public void reset(Integer v, Partition p) {
		Vertex vertex = p.vertices.get(v);
		vertex.checked++;
		if (vertex.checked >= Math.max(checked, p.edges.get(v).size())) {
			vertex.dis = 10000;
		} else {
			vertex.notInGroup = p.edges.get(v).size();

			/*
			 * for (int neighbor : p.edges.get(v)) { if(p.vertices.get(neighbor).dis!=10000)
			 * { vertex.notInGroup++; } }
			 */
		}
		// if (p.dis.get(v) != 0)
		// p.dis.put(v, 0d/*1.0 / (p.vertices.get(v).dis + 1)*/);
	}

	class repartitionThread implements Runnable {

		private Partition p;

		public repartitionThread(Partition p) {
			// TODO Auto-generated constructor stub
			this.p = p;
		}

		@Override
		public void run() {
			// System.out.println(p.id + ";" + count+"start");
			// TODO Auto-generated method stub
			// super.run();
			// while (count == 0 || Tr < 1) {
			// System.out.println(this.getName()+" start finding Groups");
			// for (Vertex v : p.subgraph().getVertices().values()) {
			/*
			 * for (Iterator<Vertex> it = p.vertices.values().iterator(); it.hasNext();) {
			 * Vertex v = it.next(); if (v.replicaSet.size() > 1 && hasLock(v, p)) {
			 * vertexCandidates.add(v); } }
			 */
			Queue<Group> groupList = new PriorityQueue<Group>(new Comparator<Group>() {

				@Override
				public int compare(Group o1, Group o2) {
					// TODO Auto-generated method stub
					if (o1.gain > o2.gain) {
						return -1;
					} else if (o1.gain < o2.gain) {
						return 1;
					}
					return 0;
				}
			});
			List<Vertex> vertexCandidates = markDis(p);
			// long t1=System.currentTimeMillis();
			for (Vertex vertex : vertexCandidates) {
				if (vertex.dis != 10000) {
					Group group = findGroup(vertex, p);

					if (group != null) {
						groupList.add(group);
					}
				}
			}
			for (Vertex v : p.vertices.values()) {
				v.dis = 10000;
				v.checked = 0;
			}
			AllgroupList.put(p.id, groupList);
			/*
			 * for (Group group : groupList) { for (Edge e : group.edgesToSend) {
			 * p.addEdgeToPartition(e); } }
			 */
			// System.out.println("partition"+p.id+" finish iterator "+count);
			// long t2=System.currentTimeMillis();
			// System.out.println("partition"+p.id+": "+(t2-t1));
			try {
				cb.await();
			} catch (InterruptedException | BrokenBarrierException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// System.out.println(this.getName()+" finish migration Groups");
			/*
			 * try { cb.await(); } catch (BrokenBarrierException e) { // TODO Auto-generated
			 * catch block e.printStackTrace(); } catch (InterruptedException e) { // TODO
			 * Auto-generated catch block e.printStackTrace(); }
			 */
			// System.out.println(p.id + ";" + count+"finish");
		}
	}

	public static void main(String[] args) throws IOException, ExecutionException {
		// TODO Auto-generated method stub
		String[] datas = { /*"dblp", "roadNet-PA", "youtube",*/ "lj", "orkut" };
		int[] ps = { 8/* , 4, 8, 16 */ };
		String[] trackers = { "hdrf"/* , "random", "Ne" */ };
		for (String string : datas) {
			String data = string;
			for (int i : ps) {
				int pNum = i;
				for (String string2 : trackers) {
					String tracker = string2;
					float balance = 1.1f;// 平衡系数
					int edgeChecked = 4;
					String inputOrder = "RDM";// 输入顺序，随机或BFS
					int mgz = 20;// 组团边最大边数量
					float gainRF = 0.0005f;
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
					if (tracker == "Ne") {
						iEP = "D:\\dataset\\" + tracker + "\\" + data + "2-partition-" + pNum;
					} else {
						iEP = "D:\\dataset\\" + tracker + "\\" + data + "2-" + inputOrder + "-partition-" + pNum;
					}
					System.out.println("-----------------" + data + "-" + pNum + "-" + tracker + "-----------------");
					GAER partitioner = new GAER(pNum, balance, vNum, iEP, mgz, edgeChecked, gainRF);
					partitioner.load();
					try {
						partitioner.repartitioning();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

				}
			}
		}
	}

}
