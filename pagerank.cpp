#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <boost/algorithm/string.hpp>
#include <pthread.h>
#include <atomic>

using namespace std;

// index for the pagerank labels,
// note that every occurrence of "CURRENT" and "NEXT" will be replaced with 0 and 1
// you can disable this and simply use 0 and 1 if you find it easier to name your variable this way
#define CURRENT 0
#define NEXT 1

pthread_barrier_t barrier;
struct timespec start, finish;
int num_threads = 4;

class CsrGraph {
// You shouldn't need to modify this class except for adding essential locks and relevant methods,
// but if you find it useful for your code you can add some minor helper functions.
// Make sure not to store too much additional information here.
// The whole point of using CSR is that it utilizes the space efficiently.

private:
  // use an array of struct for labels
  struct Node {
    // labels[0] is the current label, labels[1] is the next label
    double labels[2];
    atomic<double> labels_atomic[2];
  };

  // helper struct for COO
  struct Edge {
    int src, dst, weight;
    Edge(int s, int d, int w): src(s), dst(d), weight(w) { }
  };

  // number of nodes and edges in the graph
  int num_nodes;
  int num_edges;

  // csr representation of the graph (same as in the slides)
  int* rp;
  int* ci;
  int* ai;

  // Fine-grained locks for each node
  pthread_mutex_t* node_locks;
  pthread_spinlock_t* node_spinlocks;

  // pagerank labels
  Node* node;

  bool is_allocated;

public:
  CsrGraph() {
    // no dynamic allocated space used
    is_allocated = false;
  }

  CsrGraph(ifstream& f_dimacs) {
    // parse the input file to COO format and set number of nodes and edges
    vector<Edge> edges = parse_dimacs(f_dimacs);

    // convert COO to CSR and store the representation in the graph
    gen_csr(edges);

    // Initialize locks
    node_locks = new pthread_mutex_t [num_nodes + 1];
    node_spinlocks = new pthread_spinlock_t [num_nodes + 1];
    for (int i = 0; i <= num_nodes; i++) {
      pthread_mutex_init(&node_locks[i], NULL);
      pthread_spin_init(&node_spinlocks[i], 0);
      node[i].labels_atomic[CURRENT].store(0.0);
      node[i].labels_atomic[NEXT].store(0.0);
    }

    // remind you to delete the space
    is_allocated = true;
  }

  ~CsrGraph() {
    if (!is_allocated) {
      return;
    }

    // clean up dynamic arrays
    delete [] node;
    delete [] rp;
    delete [] ci;
    delete [] ai;
  }

  vector<Edge> parse_dimacs(ifstream& f_dimacs) {
    string line;
    vector<Edge> edges;

    // parse dimacs
    while (getline(f_dimacs, line)) {
      // skip comments
      if ('c' == line[0]) {
        continue;
      }

      // split the line into tokens
      vector<string> token;
      boost::split(token, line, [] (char c) { return c == ' '; });

      // p <problem_type> <num_nodes> <num_edges>
      // get the number of nodes in the graph
      if ("p" == token[0]) {
        num_nodes = stoi(token[2]);
        num_edges = stoi(token[3]);
      }

      // a <src> <dst> <weight>
      else if ("a" == token[0]) {
        edges.push_back(Edge(stoi(token[1]), stoi(token[2]), stoi(token[3])));
      }
    }

    return edges;
  }

  void gen_csr(vector<Edge>& edges) {
    // sort the edges
    sort(edges.begin(), edges.end(), 
        [] (const Edge& e1, const Edge& e2) { 
          return (e1.src < e2.src) ? true :
                (e1.src == e2.src) ? (e1.dst < e2.dst) : false;
        });

    // handle duplicate in the edges
    vector<Edge> good_edges;
    int edge_size = edges.size();
    int src = edges[0].src;
    int dst = edges[0].dst;
    int max_weight = edges[0].weight;
    for (int i = 1; i < edge_size; i++) {
      if (edges[i].src == edges[i - 1].src && edges[i].dst == edges[i - 1].dst) {
        // duplicate edge
        max_weight = edges[i].weight > max_weight ? edges[i].weight : max_weight;
      } 
      else {
        // add the edge with max_weight to good_edges
        good_edges.push_back (Edge(edges[i - 1].src, edges[i - 1].dst, max_weight));

        // initialize max_weight
        max_weight = edges[i].weight;
      }
    }
    // deal with the last edge
    good_edges.push_back(Edge(edges[edge_size - 1].src, edges[edge_size - 1].dst, max_weight));

    // get the actual number of edges
    num_edges = good_edges.size();

    // allocate space for dynamic arrays
    node = new Node [num_nodes + 2];
    rp = new int [num_nodes + 2];
    ci = new int [num_edges + 2];
    ai = new int [num_edges + 2];

    // rp[0] = 1 to be consistent with rp[1] = 1 (This affects the number of edges we get for node 1)
    rp[0] = 1;
    ci[0] = 0;
    ai[0] = 0;
    good_edges.insert(good_edges.begin(), Edge(rp[0], ci[0], ai[0]));
    
    // loop through each edge in COO and construct CSR representation
    int cur_node = 0;
    int cur_edge = 1;
    while (cur_edge <= num_edges) {
      const Edge& e = good_edges[cur_edge];
      // update cur_node and its coresponding rp when current source node is larger than current node
      if (e.src > cur_node) {
        cur_node++;
        rp[cur_node] = cur_edge;
      }
      // copy the destination and weight to ci and ai otherwise
      else {
        ci[cur_edge] = e.dst;
        ai[cur_edge] = e.weight;
        cur_edge++;
      }
    }

    // set the remaining nodes without edges
    for (int i = cur_node + 1; i <= num_nodes + 1; i++) {
      rp[i] = num_edges + 1;
    }
    ci[num_edges + 1] = num_nodes + 1;
    ai[num_edges + 1] = 0;
  }

  int node_size() {
    // return the number of nodes
    return num_nodes; 
  }
  int size_edges() {
    // return the number of edges
    return num_edges; 
  }

  int edge_begin(int n) {
    // return the starting index for out-going edges from node n
    return rp[n];
  }

  int edge_end(int n) {
    // return the ending index for out-going edges from node n
    return rp[n + 1]; 
  }

  double get_label(int n, int which) {
    // return the current or next label for node n
    return node[n].labels[which]; 
  }

  double get_label_atomic(int n, int which) {
    // return the current or next label for node n
    return node[n].labels_atomic[which].load(); 
  }

  void set_label(int n, int which, double label) {
    // set the current or next label for node n
    node[n].labels[which] = label;
  }

  void set_label_atomic(int n, int which, double label) {
    double cur = node[n].labels_atomic[which];
    while (!node[n].labels_atomic[which].compare_exchange_weak(cur, label));
  }

  void increment_label_atomic(int n, int which, double label) {
    double t = node[n].labels_atomic[which];
    while (!node[n].labels_atomic[which].compare_exchange_weak(t, t + label));
  }

  int get_out_degree(int n) {
    // return the number of out-degrees for node n
    return rp[n + 1] - rp[n]; 
  }

  int get_edge_dst(int e) {
    // return the destination node in edge e (index for ci)
    return ci[e];
  }

  int get_edge_weight(int e) {
    // return the weight of edge e (index for ai)
    return ai[e];
  }

  void acquire_lock(int n) {
    // acquire the lock for node n
    pthread_mutex_lock(&node_locks[n]);
  }

  void release_lock(int n) {
    // release the lock for node n
    pthread_mutex_unlock(&node_locks[n]);
  }

  void acquire_spinlock(int n) {
    // acquire the spinlock for node n
    pthread_spin_lock(&node_spinlocks[n]);
  }

  void release_spinlock(int n) {
    // release the spinlock for node n
    pthread_spin_unlock(&node_spinlocks[n]);
  }
};

void reset_next_label(CsrGraph* g, const double damping) {
  // Modify this function in any way you want to make pagerank parallel
  int num_nodes = g->node_size();
  for (int n = 1; n <= num_nodes; n++) {
    g->set_label(n, NEXT, (1.0 - damping) / (double) num_nodes);
  }
}

/**
 * Resets the next label for each node in parallel.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The CsrGraph object representing the graph.
 * @param damping The damping factor for the PageRank algorithm.
 */
void reset_next_label_parallel(int thread_id, CsrGraph* g, const double damping) {
  int num_nodes = g->node_size();

  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    g->set_label(n, NEXT, (1.0 - damping) / (double) num_nodes);
  }
}

/**
 * Resets the next label of nodes in a graph using atomic operations.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The graph to operate on.
 * @param damping The damping factor for the PageRank algorithm.
 */
void reset_next_label_atomic(int thread_id, CsrGraph* g, const double damping) { 
  int num_nodes = g->node_size();

  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    g->set_label_atomic(n, NEXT, (1.0 - damping) / (double) num_nodes);
  }
}

int iteration = 0;

bool is_converged(CsrGraph* g, const double threshold) {
  // Modify this function in any way you want to make pagerank parallel
  for (int n = 1; n <= g->node_size(); n++) {
    const double cur_label = g->get_label(n, CURRENT);
    const double next_label = g->get_label(n, NEXT);
    const double percent_change = fabs(1 - (next_label / cur_label));
    if (percent_change > threshold) {
      return false;
    }
  }  
  return true;
}

/**
 * Checks if the PageRank algorithm has converged for a specific thread in parallel execution.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The CsrGraph object representing the graph.
 * @param threshold The convergence threshold.
 * @return True if the algorithm has converged, False otherwise.
 */
bool is_converged_parallel(int thread_id, CsrGraph* g, const double threshold) {
  // if (thread_id == 0 && ++iteration == 27) {
  //   cout << "Should stop here" << endl;
  // }
  // Modify this function in any way you want to make pagerank parallel
  int num_nodes = g->node_size();
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    const double cur_label = g->get_label(n, CURRENT);
    const double next_label = g->get_label(n, NEXT);
    const double percent_change = fabs(1 - (cur_label / next_label));
    if (percent_change > threshold) {
      return false;
    }
  }  
  return true;
}

/**
 * Checks if the PageRank algorithm has converged for a specific thread.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The CsrGraph object representing the graph.
 * @param threshold The convergence threshold.
 * @return True if the algorithm has converged, false otherwise.
 */
bool is_converged_atomic(int thread_id, CsrGraph* g, const double threshold) {
  // Modify this function in any way you want to make pagerank parallel
  int num_nodes = g->node_size();
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    const double cur_label = g->get_label_atomic(n, CURRENT);
    const double next_label = g->get_label_atomic(n, NEXT);
    const double percent_change = fabs(1 - (cur_label / next_label));
    if (percent_change > threshold) {
      return false;
    }
  }  
  return true;
}

void update_current_label(CsrGraph* g) {
  // Modify this function in any way you want to make pagerank parallel
  for (int n = 1; n <= g->node_size(); n++) {
    g->set_label(n, CURRENT, g->get_label(n, NEXT));
  }
}

/**
 * Updates the current label of nodes in parallel.
 * 
 * This function updates the current label of nodes in a CsrGraph object in parallel.
 * Each thread is assigned a specific range of nodes to update.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The CsrGraph object containing the nodes.
 */
void update_current_label_parallel(int thread_id, CsrGraph* g) {
  int num_nodes = g->node_size();
  // Modify this function in any way you want to make pagerank parallel
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    g->set_label(n, CURRENT, g->get_label(n, NEXT));
  }
}

/**
 * Updates the current label of nodes in a parallel manner.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The CsrGraph object representing the graph.
 */
void update_current_label_atomic(int thread_id, CsrGraph* g) {
  int num_nodes = g->node_size();
  // Modify this function in any way you want to make pagerank parallel
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    g->set_label_atomic(n, CURRENT, g->get_label_atomic(n, NEXT));
  }
}

void scale(CsrGraph* g) {
  // Modify this function in any way you want to make pagerank parallel
  double sum = 0.0;
  for (int n = 1; n <= g->node_size(); n++) {
    sum += g->get_label(n, CURRENT);
  }
  for (int n = 1; n <= g->node_size(); n++) {
    g->set_label(n, CURRENT, g->get_label(n, CURRENT) / sum);
  }
}

double o_sum;
pthread_mutex_t o_sum_mutex;

/**
 * Scales the labels of the nodes in a CsrGraph object in a parallel manner.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The CsrGraph object containing the nodes and their labels.
 */
void scale_mutex(int thread_id, CsrGraph* g) {
  // Modify this function in any way you want to make pagerank parallel
  double sum = 0.0;
  int num_nodes = g->node_size();
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    sum += g->get_label(n, CURRENT);
  }
  pthread_mutex_lock(&o_sum_mutex);
  o_sum += sum;
  pthread_mutex_unlock(&o_sum_mutex);

  pthread_barrier_wait(&barrier);

  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    g->set_label(n, CURRENT, g->get_label(n, CURRENT) / o_sum);
  }
}

pthread_spinlock_t o_sum_spinlock;

/**
 * Scales the labels of the nodes in a graph using a spinlock for parallelization.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The CsrGraph object representing the graph.
 */
void scale_spinlock(int thread_id, CsrGraph* g) {
  // Modify this function in any way you want to make pagerank parallel
  double sum = 0.0;
  int num_nodes = g->node_size();
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    sum += g->get_label(n, CURRENT);
  }
  pthread_spin_lock(&o_sum_spinlock);
  o_sum += sum;
  pthread_spin_unlock(&o_sum_spinlock);

  pthread_barrier_wait(&barrier);

  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    g->set_label(n, CURRENT, g->get_label(n, CURRENT) / o_sum);
  }
}

atomic<double> o_sum_atomic{0.0};

/**
 * Scales the labels of the nodes in a CsrGraph object in a parallel manner.
 * 
 * @param thread_id The ID of the current thread.
 * @param g The CsrGraph object containing the nodes and labels.
 */
void scale_atomic(int thread_id, CsrGraph* g) {
  // Modify this function in any way you want to make pagerank parallel
  double sum = 0.0;
  int num_nodes = g->node_size();
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    sum += g->get_label_atomic(n, CURRENT);
  }

  double o_sum = o_sum_atomic;
  while (!o_sum_atomic.compare_exchange_weak(o_sum, o_sum + sum));

  pthread_barrier_wait(&barrier);

  double f_sum = o_sum_atomic;
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    g->set_label_atomic(n, CURRENT, g->get_label_atomic(n, CURRENT) / f_sum);
  }
}

struct thread_data {
  int thread_id;
  CsrGraph* g;
  double threshold;
  double damping;
};

bool converged = true;
pthread_spinlock_t converged_spinlock;

void init_spin() {
  pthread_spin_init(&converged_spinlock, 0);
  pthread_spin_init(&o_sum_spinlock, 0);
}

void* compute_pagerank_spinlock(void* threadarg) {
  // You have to divide the work and assign it to threads to make this function parallel
  struct thread_data* my_data = (struct thread_data*) threadarg;
  CsrGraph* g = my_data->g;
  double threshold = my_data->threshold;
  double damping = my_data->damping;
  int thread_id = my_data->thread_id;

  // initialize
  bool convergence = false;
  int num_nodes = g->node_size();
  for (int n = thread_id; n <= num_nodes; n+=num_threads) {
    g->set_label(n, CURRENT, 1.0 / num_nodes);
  }

  pthread_barrier_wait(&barrier);

  if (thread_id == 0) 
    clock_gettime(CLOCK_MONOTONIC, &start);

  do {
    pthread_barrier_wait(&barrier);

    if (thread_id == 0) 
      converged = true;
    // reset next labels
    reset_next_label_parallel(thread_id, g, damping);

    pthread_barrier_wait(&barrier);

    // apply current node contribution to others
    for (int n = thread_id; n <= num_nodes; n+=num_threads) {
      // compute total out weight
      double out_weight = 0.0;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        out_weight += g->get_edge_weight(e);
      }
      out_weight = out_weight > 0.0 ? out_weight : 1.0;

      double my_contribution = damping * g->get_label(n, CURRENT) / out_weight;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        int dst = g->get_edge_dst(e);
        g->acquire_spinlock(dst);
        g->set_label(dst, NEXT, g->get_label(dst, NEXT) + my_contribution);
        g->release_spinlock(dst);
      }
    }

    pthread_barrier_wait(&barrier);

    // check the change across successive iterations to determine convergence
    bool my_converged = is_converged_parallel(thread_id, g, threshold);
    pthread_spin_lock(&converged_spinlock);
    converged = converged && my_converged;
    pthread_spin_unlock(&converged_spinlock);

    // update current labels
    update_current_label_parallel(thread_id, g);

    pthread_barrier_wait(&barrier);
  } while(!converged);

  if (thread_id == 0) 
    clock_gettime(CLOCK_MONOTONIC, &finish);

  // scale the sum to 1
  scale_spinlock(thread_id, g);

  return NULL;
}

pthread_mutex_t converged_mutex;

void init_mutex() {
  pthread_mutex_init(&converged_mutex, NULL);
  pthread_mutex_init(&o_sum_mutex, NULL);
}

void* compute_pagerank_mutex(void* threadarg) {
  // You have to divide the work and assign it to threads to make this function parallel
  struct thread_data* my_data = (struct thread_data*) threadarg;
  CsrGraph* g = my_data->g;
  double threshold = my_data->threshold;
  double damping = my_data->damping;
  int thread_id = my_data->thread_id;

  // initialize
  int num_nodes = g->node_size();
  for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
    g->set_label(n, CURRENT, 1.0 / num_nodes);
  }

  pthread_barrier_wait(&barrier);

  if (thread_id == 0) 
    clock_gettime(CLOCK_MONOTONIC, &start);

  do {
    pthread_barrier_wait(&barrier);

    if (thread_id == 0) 
      converged = true;

    // reset next labels
    reset_next_label_parallel(thread_id, g, damping);

    pthread_barrier_wait(&barrier);

    // apply current node contribution to others
    for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
      // compute total out weight
      double out_weight = 0.0;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        out_weight += g->get_edge_weight(e);
      }
      out_weight = out_weight > 0.0 ? out_weight : 1.0;

      double my_contribution = damping * g->get_label(n, CURRENT) / out_weight;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        int dst = g->get_edge_dst(e);
        g->acquire_lock(dst);
        g->set_label(dst, NEXT, g->get_label(dst, NEXT) + my_contribution);
        g->release_lock(dst);
      }
    }

    pthread_barrier_wait(&barrier);

    // check the change across successive iterations to determine convergence
    bool my_converged = is_converged_parallel(thread_id, g, threshold);
    pthread_mutex_lock(&converged_mutex);
    converged = converged && my_converged;
    pthread_mutex_unlock(&converged_mutex);

    // update current labels
    update_current_label_parallel(thread_id, g);

    pthread_barrier_wait(&barrier);
  } while(!converged);

  if (thread_id == 0) 
    clock_gettime(CLOCK_MONOTONIC, &finish);

  // scale the sum to 1
  scale_mutex(thread_id, g);

  return NULL;
}

atomic<bool> converged_atomic;

void init_atomic() {
  o_sum_atomic.store(0.0);
  converged_atomic.store(true);
}

void* compute_pagerank_atomic(void* threadarg) {
  // You have to divide the work and assign it to threads to make this function parallel
  struct thread_data* my_data = (struct thread_data*) threadarg;
  CsrGraph* g = my_data->g;
  double threshold = my_data->threshold;
  double damping = my_data->damping;
  int thread_id = my_data->thread_id;

  // initialize
  int num_nodes = g->node_size();
  for (int n = thread_id; n <= num_nodes; n+=num_threads) {
    g->set_label_atomic(n, CURRENT, 1.0 / num_nodes);
  }

  pthread_barrier_wait(&barrier);

  if (thread_id == 0) 
    clock_gettime(CLOCK_MONOTONIC, &start);

  do {
    pthread_barrier_wait(&barrier);

    if (thread_id == 0) 
      converged_atomic.store(true);
    // reset next labels
    reset_next_label_atomic(thread_id, g, damping);

    pthread_barrier_wait(&barrier);

    // apply current node contribution to others
    for (int n = thread_id + 1; n <= num_nodes; n+=num_threads) {
      // compute total out weight
      double out_weight = 0.0;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        out_weight += g->get_edge_weight(e);
      }
      out_weight = out_weight > 0.0 ? out_weight : 1.0;

      double my_contribution = damping * g->get_label_atomic(n, CURRENT) / out_weight;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        int dst = g->get_edge_dst(e);
        g->increment_label_atomic(dst, NEXT, my_contribution);
      }
    }

    pthread_barrier_wait(&barrier);

    // check the change across successive iterations to determine convergence
    bool converge = converged_atomic;
    while (!converged_atomic.compare_exchange_weak(converge, converge && is_converged_atomic(thread_id, g, threshold)));

    // update current labels
    update_current_label_atomic(thread_id, g);

    pthread_barrier_wait(&barrier);
  } while(!converged_atomic.load());

  if (thread_id == 0) 
    clock_gettime(CLOCK_MONOTONIC, &finish);

  // scale the sum to 1
  scale_atomic(thread_id, g);

  return NULL;
}

/**
 * Searches for the vertex range in a given graph based on the specified edge number.
 * 
 * @param g The CsrGraph object representing the graph.
 * @param start The starting index of the vertex range to search.
 * @param end The ending index of the vertex range to search.
 * @param edge_num The edge number to search for within the vertex range.
 * @return The index of the vertex range that contains the specified edge number, or the ending index if not found.
 */
int searchForVertexRange(CsrGraph* g, int start, int end, int edge_num) {
  int left = start;
  int right = end;
  int mid = (left + right) / 2;
  while (left < right) {
    if (g->edge_begin(mid) <= edge_num && g->edge_end(mid) > edge_num) {
      return mid;
    }
    else if (g->edge_end(mid) <= edge_num) {
      left = mid + 1;
    }
    else {
      right = mid;
    }
    mid = (left + right) / 2;
  }
  return end;
}

void* compute_pagerank_edge_atomic(void* threadarg) {
  // You have to divide the work and assign it to threads to make this function parallel
  struct thread_data* my_data = (struct thread_data*) threadarg;
  CsrGraph* g = my_data->g;
  double threshold = my_data->threshold;
  double damping = my_data->damping;
  int thread_id = my_data->thread_id;

  // initialize
  int num_nodes = g->node_size();
  int num_edges = g->size_edges();
  for (int n = thread_id; n <= num_nodes; n+=num_threads) {
    g->set_label_atomic(n, CURRENT, 1.0 / num_nodes);
  }

  int start_edge = thread_id == 0 ? 1 : 
                    (num_edges * (thread_id / (double) num_threads));
  int end_edge = thread_id == num_threads - 1 ? g->size_edges() + 1 : 
                    (num_edges * ((thread_id + 1) / (double) num_threads));

  int start_vertex = searchForVertexRange(g, 1, num_nodes, start_edge);
  int end_vertex = searchForVertexRange(g, 1, num_nodes, end_edge);

  pthread_barrier_wait(&barrier);

  if (thread_id == 0) 
    clock_gettime(CLOCK_MONOTONIC, &start);

  do {
    pthread_barrier_wait(&barrier);

    if (thread_id == 0) 
      converged_atomic.store(true);

    // reset next labels
    reset_next_label_atomic(thread_id, g, damping);

    pthread_barrier_wait(&barrier);

    // apply current node contribution to others
    for (int n = start_vertex; n <= end_vertex; n++) {
      // compute total out weight
      double out_weight = 0.0;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        out_weight += g->get_edge_weight(e);
      }
      out_weight = out_weight > 0.0 ? out_weight : 1.0;

      double my_contribution = damping * g->get_label_atomic(n, CURRENT) / out_weight;
      int e = n == start_vertex ? start_edge : g->edge_begin(n);
      int stop = n == end_vertex ? end_edge : g->edge_end(n);
      for (; e < stop; e++) {
        int dst = g->get_edge_dst(e);
        g->increment_label_atomic(dst, NEXT, my_contribution);
      }
    }

    pthread_barrier_wait(&barrier);

    // check the change across successive iterations to determine convergence
    bool converge = converged_atomic;
    while (!converged_atomic.compare_exchange_weak(converge, converge && is_converged_atomic(thread_id, g, threshold)));

    // update current labels
    update_current_label_atomic(thread_id, g);

    pthread_barrier_wait(&barrier);
  } while(!converged_atomic.load());

  if (thread_id == 0) 
    clock_gettime(CLOCK_MONOTONIC, &finish);

  // scale the sum to 1
  scale_atomic(thread_id, g);

  return NULL;
}

void compute_pagerank_serial(CsrGraph* g, const double threshold, const double damping) {
  // You have to divide the work and assign it to threads to make this function parallel

  // initialize
  bool convergence = false;
  int num_nodes = g->node_size();
  for (int n = 1; n <= num_nodes; n++) {
    g->set_label(n, CURRENT, 1.0 / num_nodes);
  }

  clock_gettime(CLOCK_MONOTONIC, &start);

  do {
    // reset next labels
    reset_next_label(g, damping);

    // apply current node contribution to others
    for (int n = 1; n <= num_nodes; n++) {
      // compute total out edge weight
      double out_weight = 0.0;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        out_weight += g->get_edge_weight(e);
      }
      out_weight = out_weight > 0.0 ? out_weight : 1.0;

      double my_contribution = damping * g->get_label(n, CURRENT) / out_weight;
      for (int e = g->edge_begin(n); e < g->edge_end(n); e++) {
        int dst = g->get_edge_dst(e);
        g->set_label(dst, NEXT, g->get_label(dst, NEXT) + my_contribution);
      }
    }

    // check the change across successive iterations to determine convergence
    convergence = is_converged(g, threshold);

    // update current labels
    update_current_label(g);
  } while(!convergence);

  clock_gettime(CLOCK_MONOTONIC, &finish);

  // scale the sum to 1
  scale(g);
}

void compute_pagerank(int type, CsrGraph* g, const double threshold, const double damping) {
  if (type == 0) {
    compute_pagerank_serial(g, threshold, damping);
    return;
  }

  // Create threads
  pthread_t threads[num_threads];
  int thread_ids[num_threads];

  // Create barrier
  pthread_barrier_init(&barrier, NULL, num_threads);

  // Create function pointer to function with desired type of synchronization
  void* (*compute_pagerank_func)(void*);
  switch(type) {
    case 1:
      init_mutex();
      compute_pagerank_func = compute_pagerank_mutex;
      break;
    case 2:
      init_spin();
      compute_pagerank_func = compute_pagerank_spinlock;
      break;
    case 3:
      init_atomic();
      compute_pagerank_func = compute_pagerank_atomic;
      break;
    case 4:
      init_atomic();
      compute_pagerank_func = compute_pagerank_edge_atomic;
      break;
    default:
      init_mutex();
      compute_pagerank_func = compute_pagerank_mutex;
      break;
  }

  for (int i = 0; i < num_threads; i++) {
    thread_ids[i] = i;
    struct thread_data* data = new thread_data {i, g, threshold, damping};
    pthread_create(&threads[i], NULL, compute_pagerank_func, (void*) data);
  }

  for (int i = 0; i < num_threads; i++) {
    pthread_join(threads[i], NULL);
  }

  if (type == 3 || type == 4) {
    for (int n = 1; n <= g->node_size(); n++) {
      g->set_label(n, CURRENT, g->get_label_atomic(n, CURRENT));
    }
  }
}

void sort_and_print_label(CsrGraph* g, string out_file) {
  // You shouldn't need to change this.

  // prepare the label to be sorted
  vector<pair<int, double>> label;
  for (int n = 1; n <= g->node_size(); n++) {
    label.push_back(make_pair(n, g->get_label(n, CURRENT)));
  }

  // sort the labels in descending order then node number in ascending order
  sort(label.begin(), label.end(), 
      [] (const pair<int, double>& v1, const pair<int, double>& v2) {
        return (v1.second > v2.second) ? true :
                (v1.second == v2.second) ? (v1.first < v2.first) : false;
      });

  // Create output stream
  ofstream out_stream(out_file);

  // Compute time in nanoseconds
  long total_time = (finish.tv_sec - start.tv_sec) * 1000000000 + (finish.tv_nsec - start.tv_nsec);
  out_stream << "Time: " << total_time << " ns" << endl;

  // Print labels
  for (const pair<int, double>& v: label) {
    out_stream << v.first << " " << fixed << setprecision(6) << v.second << endl;
  }
}

int main(int argc, char *argv[]) {
  // Ex: ./pagerank road-NY.dimacs road-NY.txt
  if (argc < 5) {
    cerr << "Usage: " << argv[0] << " <input.dimacs> <output_filename>\n";
    return 0;
  }

  // make sure the input argument is valid
  ifstream f_dimacs(argv[1]);
  if (!f_dimacs) {
    cerr << "Failed to open " << argv[1] << endl;
    return 0;
  }

  // construct the CSR graph
  CsrGraph* g = new CsrGraph(f_dimacs);

  // define important constants for pagerank
  num_threads = atoi(argv[4]);
  int type = atoi(argv[3]);
  const double threshold = 0.01;
  const double damping = 0.85;

  // compute the pagerank using push-style method
  compute_pagerank(type, g, threshold, damping);

  // sort and print the labels to the output file
  sort_and_print_label(g, argv[2]);

  // delete the allocated space to the graph avoid memory leak
  delete g;

  return 0;
}