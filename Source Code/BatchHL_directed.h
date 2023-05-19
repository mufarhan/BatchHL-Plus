#ifndef DHGHWAY_LABELING_H_
#define DHGHWAY_LABELING_H_

#include <sys/time.h>
#include <cmath>
#include <iostream>
#include <thread>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <fstream>

#include "two_layer_queue.h"

//
// NOTE: Unweighted and directed graphs are supported.
//

class DHighwayLabelling {
public:
  // Constructs an index from a graph, given as a list of edges.
  DHighwayLabelling(std::string filename, int k);
  DHighwayLabelling();
  ~DHighwayLabelling();

  void ConstructHighwayLabelling(int i, int topk[], int dir);
  void BuildIndex(int topk[]);

  void UpdateLabelling(std::string filename);
  void BHL_Plus(std::vector<std::pair<std::string, std::pair<int, int> > > updates[], int &sum);
  void BHL(std::vector<std::pair<std::string, std::pair<int, int> > > updates[], int &sum);
  bool prunable(int i, int u, int *temp, int *A, int dir);

  void BatchSearch_Deletion(int i, int *A, int *AA, int *AFF_VERTS, int &c, std::vector<std::pair<std::string, std::pair<int, int> > > &updates, int dir);
  void BatchSearch_Combined(int i, int *A, int *temp, int *AFF_VERTS, int &c, std::vector<std::pair<std::string, std::pair<int, int> > > &updates, int dir);
  void BatchRepair(int i, int *A, int *temp, int *AFF_VERTS, int &c, int dir);

  void deallocate();
  void PrintLabelling();

  void Mutate(std::string filename_1, std::string filename_2);

  int ldPair(int d, bool landmark);
  int ldDist(int ld);
  int ldAdd(int ld, bool landmark);
  int InfoEmpty(int d, bool landmark_flag);
  int InfoAdd(int info, bool is_landmark);
  int addTwoFlags(int d, bool landmark_flag, bool stable_flag);

  void SelectRandomPairs(std::string filename);
  void SelectUpdatesRealWorld(std::string filename);

  void SelectLandmarks_HD(int topk[]);
  int LabellingSize();

  int query(int r, int v, int dir);
  int min(int a, int b);
  int max(int a, int b);

  void RemoveLandmarks(int topk[]);
  void AddLandmarks(int topk[]);

  void experiment_1(std::string index_file, std::string batch_file, int batch_size, int topk[]);

  int HC_UB_naive(int s, int t);
  int HC_LB_naive(int s, int t);
  void QueryDistance(std::string pairs, std::string output);

  void storeLabelling(std::string filename);
  void loadLabelling_Full(std::string filename, int topk[]);
  void loadLabelling_Pruned(std::string filename);

private:
  int V;  // total number of vertices
  long E; // total number of edges
  int K; // total number of landmarks

  int **distances[2], **distances_1[2];
  int **highway[2], **highway_1[2];
  int **vertices[2];
  int *C[2];
  std::vector<std::vector<int> > adj[2];
  std::unordered_map<int, std::vector<int> > adj1;
  std::map<int, int> landmarks;

  double GetCurrentTimeSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
  }

  long GetCurrentTimeMicroSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t) 1e6 + tv.tv_usec;
  }

  long GetCurrentTimeMilliSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000LL + tv.tv_usec / 1000;
  }

  // Statistics
  double time_, time_querying_sec_;
  long time_querying_microsec_, time_querying_millisec_;
};

DHighwayLabelling::DHighwayLabelling() { }

DHighwayLabelling::~DHighwayLabelling() { }

void DHighwayLabelling::deallocate() {

  for(int dir = 0; dir < 2; dir++) {
    for(int i = 0; i < V; i++) {
      delete [] distances[dir][i];
      delete [] distances_1[dir][i];
    }
    delete [] distances[dir];
    delete [] distances_1[dir];

    for(int i = 0; i < K; i++) {
      delete [] highway[dir][i];
      delete [] highway_1[dir][i];
    }
    delete [] highway[dir];
    delete [] highway_1[dir];
  }
}

DHighwayLabelling::DHighwayLabelling(std::string filename, int k) {
  V = 0; E = 0; K = k;

  std::ifstream ifs(filename);
  if (ifs.is_open()) {

    std::vector<std::pair<int, int> > es;
    std::unordered_map<int, int> vertex2id;
    for (int u, v; ifs >> u >> v;) {
      if (vertex2id.count(u) == 0) vertex2id[u] = V++;
      if (vertex2id.count(v) == 0) vertex2id[v] = V++;
      u = vertex2id[u];
      v = vertex2id[v];
      if (u != v)
        es.emplace_back(u, v);
    }

    std::sort(es.begin(), es.end());
    es.erase(std::unique(es.begin(), es.end()), es.end());
    ifs.close();

    for (int dir = 0; dir < 2; dir++)
      adj[dir].resize(V);

    for (const auto &e : es) {
      if (e.first == e.second) continue; 
      adj[0][e.first].push_back(e.second);
      adj[1][e.second].push_back(e.first);
    }

    for (int dir = 0; dir < 2; dir++) {
      for (int v = 0; v < V; v++) {
        std::sort(adj[dir][v].begin(), adj[dir][v].end());
        adj[dir][v].erase(std::unique(adj[dir][v].begin(), adj[dir][v].end()), adj[dir][v].end());
        adj[dir][v].shrink_to_fit();
      }
    }

    for (int dir = 0; dir < 2; dir++) {
      for(int v = 0; v < V; v++)
        E += adj[dir][v].size();
    }
    std::cout << "V : " << V << " E : " << E << std::endl << std::endl;
  
  } else
    std::cout << "Unable to open file" << std::endl;
}

void DHighwayLabelling::PrintLabelling() {

  for(int dir = 0; dir < 2; dir++) {
    std::cout << std::endl << "Direction " << dir << std::endl;
    for (int i = 0; i < V; i++) {
      std::cout << i << " -> ";
      for (int j = 0; j < K; j++) {
        if(distances_1[dir][i][j] != 999)
          std::cout<< "(" << j << ", " << (int) distances_1[dir][i][j] << ") ";
      }
      std::cout<<std::endl;
    }
  }

  for(int dir = 0; dir < 2; dir++) {
    std::cout << std::endl << "Direction " << dir << std::endl;
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < K; j++)
        std::cout<< (int) highway_1[dir][i][j] << " ";
      std::cout<<std::endl;
    }
  }
}

int DHighwayLabelling::LabellingSize() {
  long size = 0;

  for(int dir = 0; dir < 2; dir++) {
    for (int i = 0; i < V; i++) {
      for (int j = 0; j < K; j++) {
        if(distances_1[dir][i][j] != 999)
          size++;
      }
    }
  }
  
  return (V + 2 *size) / (1024 * 1024);
}

void DHighwayLabelling::ConstructHighwayLabelling(int i, int topk[], int dir) {

  int *P = new int[V];
  for(int j = 0; j < V; j++)
    P[j] = 999;

  std::queue<int> que[2];
  que[0].push(topk[i]); que[0].push(-1);
  distances[dir][topk[i]][i] = 0; distances_1[dir][topk[i]][i] = 0; P[topk[i]] = 0; int use = 0;
  while (!que[0].empty()) {
    int u = que[use].front();
    que[use].pop();

    if(u == -1) {
      use = 1 - use;
      que[use].push(-1);
      continue;
    }

    for (int w : adj[dir][u]) {
      if (P[w] == 999) {
        P[w] = P[u] + 1;
        if(use == 1 || landmarks.count(w) > 0)
          que[1].push(w);
        else {
          que[0].push(w);
          distances[dir][w][i] = P[w];
          distances_1[dir][w][i] = P[w];
        }
      }
    }
  }

  for(int j = 0; j < K; j++) {
    highway[dir][i][j] = P[topk[j]];
    highway_1[dir][i][j] = P[topk[j]];
  }

  delete [] P;
}

void DHighwayLabelling::BuildIndex(int topk[]) {

  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;

  /*std::ifstream ifs("../Batches_d/Indochina-B10_D.txt"); int a, b; std::string op;
  while(ifs >> op >> a >> b) {
    if(op == "EI") {
      adj[0][a].push_back(b);
      adj[1][b].push_back(a);
    } else if(op == "ED") {
      adj[0][a].erase(std::remove(adj[0][a].begin(), adj[0][a].end(), b), adj[0][a].end());
      adj[1][b].erase(std::remove(adj[1][b].begin(), adj[1][b].end(), a), adj[1][b].end());
    }
  }
  ifs.close();*/

  time_ = 0;
  for(int dir = 0; dir < 2; dir++) {
    // Initialization
    distances[dir] = new int*[V];
    distances_1[dir] = new int*[V];
    for(int i = 0; i < V; i++) {
      distances[dir][i] = new int[K];
      distances_1[dir][i] = new int[K];
      for(int j = 0; j < K; j++) {
        distances[dir][i][j] = 999;
        distances_1[dir][i][j] = 999;
      }
    }

    highway[dir] = new int*[K];
    highway_1[dir] = new int*[K];
    for(int i = 0; i < K; i++) {
      highway[dir][i] = new int[K];
      highway_1[dir][i] = new int[K];
      for(int j = 0; j < K; j++) {
        highway[dir][i][j] = 999;
        highway_1[dir][i][j] = 999;
      }
    }

    // Start computing Highway Labelling (HL)
    double temp = -GetCurrentTimeSec();
    for (int i = 0; i < K; i++)
      ConstructHighwayLabelling(i, topk, dir);
    temp += GetCurrentTimeSec();
    time_ += temp;
  }

  std::cout << "Construction Time (sec.): " << time_ << " Labelling Size: " << LabellingSize() << " MB" << std::endl;
}

void DHighwayLabelling::UpdateLabelling(std::string filename) {
  std::ifstream ifs(filename); std::string op; int a, b;
  std::vector<std::pair<std::string, std::pair<int, int> > > updates[2];

  int *A = new int[V];
  int *temp = new int[V];
  int *AFF_VERTS = new int[V]; int c;
  for(int i = 0; i < V; i++)
    A[i] = 999;

  while(ifs >> op >> a >> b) {

    if(op == "EI") {
      adj[0][a].push_back(b);
      adj[1][b].push_back(a);

      updates[0].push_back(std::make_pair(op, std::make_pair(a, b)));
      updates[1].push_back(std::make_pair(op, std::make_pair(b, a)));
    } else if(op == "ED") {
      adj[0][a].erase(std::remove(adj[0][a].begin(), adj[0][a].end(), b), adj[0][a].end());
      adj[1][b].erase(std::remove(adj[1][b].begin(), adj[1][b].end(), a), adj[1][b].end());

      updates[0].push_back(std::make_pair(op, std::make_pair(a, b)));
      updates[1].push_back(std::make_pair(op, std::make_pair(b, a)));
    }
  }
  ifs.close(); 

  time_ = -GetCurrentTimeSec();
  for(int dir = 0; dir < 2; dir++) {
    for(int i = 0; i < K; i++) {
      BatchSearch_Combined(i, A, temp, AFF_VERTS, c = 0, updates[dir], dir);
      BatchRepair(i, A, temp, AFF_VERTS, c, dir);

      for(int j = 0; j < c; j++)
        A[AFF_VERTS[j]] = 999;
    }
  }
  time_ += GetCurrentTimeSec();

  for(int dir = 0; dir < 2; dir++) {
 
    for(a = 0; a < V; a++) {
      for(b = 0; b < K; b++)
        distances_1[dir][a][b] = distances[dir][a][b];
    }

    for(a = 0; a < K; a++) {
      for(b = 0; b < K; b++) {
        highway_1[dir][a][b] = highway[dir][a][b];
      }
    }
  }

  std::cout << "Batch Update Time (sec.): " << time_ << " Updated Labelling Size: " << LabellingSize() << std::endl;

  delete [] A; delete [] temp; delete [] AFF_VERTS;
  
}

int DHighwayLabelling::ldPair(int d, bool landmark) { return landmark ? d << 1 : (d << 1) | 1; }
int DHighwayLabelling::ldDist(int ld) { return ld >> 1; }
int DHighwayLabelling::ldAdd(int ld, bool landmark) { return landmark ? (ld + 2) & ~1 : ld + 2; }

int DHighwayLabelling::addTwoFlags(int d, bool landmark_flag, bool stable_flag) {
  return landmark_flag ? d << 2 : stable_flag ? (d << 2) | 2 : (d << 2) | 3;
}

void DHighwayLabelling::BatchSearch_Combined(int i, int *A, int *temp, int *AFF_VERTS, int &c, std::vector<std::pair<std::string, std::pair<int, int> > > &updates, int dir) {

  // computing affected vertices
  std::priority_queue<std::pair<int, std::pair<int, int> >, std::vector<std::pair<int, std::pair<int, int> > >, std::greater<std::pair<int, std::pair<int, int> > > > que1;
  std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >, std::greater<std::pair<int, int> > > que2;
  std::queue<std::pair<int, int> > que; 

  int br, beta;
  for(std::pair<std::string, std::pair<int, int> > iter : updates) {

    temp[iter.second.first] = query(i, iter.second.first, dir);
    temp[iter.second.second] = query(i, iter.second.second, dir);

    if(iter.first == "EI") {
      if(temp[iter.second.first] > temp[iter.second.second]) {
        br = ldAdd(ldPair(temp[iter.second.second], distances_1[dir][iter.second.second][i] == 999), landmarks.count(iter.second.first) > 0); beta = ldPair(temp[iter.second.first], distances_1[dir][iter.second.first][i] == 999);
        if(br < beta)
          que1.push(std::make_pair(ldPair(br, true), std::make_pair(iter.second.first, iter.second.second)));
      } else if(temp[iter.second.first] < temp[iter.second.second]) {
        br = ldAdd(ldPair(temp[iter.second.first], distances_1[dir][iter.second.first][i] == 999), landmarks.count(iter.second.second) > 0); beta = ldPair(temp[iter.second.second], distances_1[dir][iter.second.second][i] == 999);
        if(br < beta)
          que1.push(std::make_pair(ldPair(br, true), std::make_pair(iter.second.second, iter.second.first)));
      }
    } else if(iter.first == "ED") {
      if(temp[iter.second.first] > temp[iter.second.second]) {
        br = ldAdd(ldPair(temp[iter.second.second], distances_1[dir][iter.second.second][i] == 999), landmarks.count(iter.second.first) > 0); beta = ldPair(temp[iter.second.first], distances_1[dir][iter.second.first][i] == 999);
        if(br == beta)
          que2.push(std::make_pair(addTwoFlags(temp[iter.second.first], false, false), iter.second.first));
      } else if(temp[iter.second.first] < temp[iter.second.second]) {
        br = ldAdd(ldPair(temp[iter.second.first], distances_1[dir][iter.second.first][i] == 999), landmarks.count(iter.second.second) > 0); beta = ldPair(temp[iter.second.second], distances_1[dir][iter.second.second][i] == 999);
        if(br == beta)
          que2.push(std::make_pair(addTwoFlags(temp[iter.second.second], false, false), iter.second.second));
      }
    }
  }

  START:
  while (!que.empty() || !que1.empty() || !que2.empty()) {

    std::pair<int, std::pair<int, int> > p1; std::pair<int, int> p2; int d_v; int v;
    if(!que.empty() && !que1.empty() && !que2.empty()) {
      if(que.front().first < que1.top().first && que.front().first < que2.top().first) {
        p2 = que.front(); que.pop();
        d_v = p2.first; v = p2.second;
      } else if(que1.top().first < que.front().first && que1.top().first < que2.top().first) {
        p1 = que1.top(); que1.pop();
        d_v = p1.first; v = p1.second.first;
        if(A[p1.second.second] != 999)
          goto START;
      } else if (que2.top().first < que.front().first && que2.top().first < que1.top().first) {
        p2 = que2.top(); que2.pop();
        d_v = p2.first; v = p2.second;
      }
    } else if(!que1.empty() && !que2.empty()) {
      if(que1.top().first < que2.top().first) {
        p1 = que1.top(); que1.pop();
        d_v = p1.first; v = p1.second.first;
        if(A[p1.second.second] != 999)
          goto START;
      } else {
        p2 = que2.top(); que2.pop();
        d_v = p2.first; v = p2.second;
      }
    } else if(!que.empty() && !que1.empty()) {
      if(que1.top().first < que.front().first) {
        p1 = que1.top(); que1.pop();
        d_v = p1.first; v = p1.second.first;
        if(A[p1.second.second] != 999)
          goto START;
      } else {
        p2 = que.front(); que.pop();
        d_v = p2.first; v = p2.second;
      }
    } else if(!que.empty() && !que2.empty()) {
      if(que.front().first < que2.top().first) {
        p2 = que.front(); que.pop();
        d_v = p2.first; v = p2.second;
      } else {
        p2 = que2.top(); que2.pop();
        d_v = p2.first; v = p2.second;
      }
    } else if(!que1.empty()) {
      p1 = que1.top(); que1.pop();
      d_v = p1.first; v = p1.second.first;
      if(A[p1.second.second] != 999)
        goto START;
    } else if(!que2.empty()) {
      p2 = que2.top(); que2.pop();
      d_v = p2.first; v = p2.second;
    } else {
      p2 = que.front(); que.pop();
      d_v = p2.first; v = p2.second;
    }

    if(A[v] == 999) {
      if (d_v >> 1 == (temp[v] << 1) + 1) {
        // check for path of equal (landmark) length
        for(int w : adj[1-dir][v]) {
          if(A[w] != 1001) {
            int ld_w;
            if(A[w] != 999)
              ld_w = A[w];
            else {
              ld_w = ldPair(query(i, w, dir), distances_1[dir][w][i] == 999);
            }

            if(d_v >> 2 == (ld_w >> 1) + 1) {
              d_v &= ~1;
              if(ldPair(temp[v], distances_1[dir][v][i] == 999) == ldAdd(ld_w, landmarks.count(v) > 0))
                goto START;
            }
          }
        }
      }
      A[v] = 1000; // ID 1000 is for tracking V_{AFF}+
      AFF_VERTS[c] = v; c++;

      if((d_v & 1) == 0) {
        A[v] = d_v >> 1;

        for(int w : adj[dir][v]) {
	  temp[w] = query(i, w, dir);
          br = ldAdd(ldPair(temp[v], distances_1[dir][v][i] == 999), landmarks.count(w) > 0); beta = ldPair(temp[w], distances_1[dir][w][i] == 999); int br1 = ldAdd(A[v], landmarks.count(w) > 0);
          if(br == beta && beta < br1 || br1 < beta)
            que.push(std::make_pair(ldPair(br1, true), w));
        }

      } else {
	A[v] = 1001; // ID 1001 is for tracking V_{UNSTABLE}
        for(int w : adj[dir][v]) {
          temp[w] = query(i, w, dir);
          br = ldAdd(ldPair(temp[v], distances_1[dir][v][i] == 999), landmarks.count(w) > 0); beta = ldPair(temp[w], distances_1[dir][w][i] == 999);
          if((d_v >> 2) + 1 < temp[w] || br == beta)
            que.push(std::make_pair(addTwoFlags((d_v >> 2) + 1, false, false), w));
        }
      }
    }
  }
}

void DHighwayLabelling::BatchRepair(int i, int *A, int *temp, int *AFF_VERTS, int &c, int dir) {

  // computing boundary vertices
  std::vector<std::pair<int, int> > V_aff;
  for(int j = 0; j < c; j++) {
    temp[AFF_VERTS[j]] = 999;
    for (int w : adj[1 - dir][AFF_VERTS[j]]) {
      if(A[w] == 999) {
	temp[w] = query(i, w, dir);
        temp[AFF_VERTS[j]] = min(temp[AFF_VERTS[j]], temp[w] + 1);
      }
    }

    if(temp[AFF_VERTS[j]] != 999)
      V_aff.push_back(std::make_pair(temp[AFF_VERTS[j]], AFF_VERTS[j]));
    landmarks.count(AFF_VERTS[j])>0?highway[dir][i][landmarks[AFF_VERTS[j]]] = 999:distances[dir][AFF_VERTS[j]][i] = 999;
  }

  // updating the labelling
  if(V_aff.size() > 0) {
    std::sort(V_aff.begin(), V_aff.end());

    std::queue<int> quee[2]; int use = 0;
    int d = V_aff[0].first; int x = 0;
    while(x < V_aff.size() && V_aff[x].first == d) {
      if(prunable(i, V_aff[x].second, temp, A, dir)){
        quee[1].push(V_aff[x].second);
      } else {
        distances[dir][V_aff[x].second][i] = temp[V_aff[x].second];
        quee[0].push(V_aff[x].second);
      }
      A[V_aff[x].second] = 999;
      x++;
    }

    quee[use].push(-1);
    while (!quee[0].empty() || x < V_aff.size()) {
      int u = quee[use].front();
      quee[use].pop();
       if(u == -1) {
        if(use == 0) { d++;
          while(x < V_aff.size() && V_aff[x].first == d) {
            if(A[V_aff[x].second] != 999) {
              if(prunable(i, V_aff[x].second, temp, A, dir)) {
                quee[1].push(V_aff[x].second);
              } else {
                distances[dir][V_aff[x].second][i] = temp[V_aff[x].second];
                quee[0].push(V_aff[x].second);
              }
              A[V_aff[x].second] = 999;
            }
            x++;
          }
        }
        use = 1 - use;
        quee[use].push(-1);
        continue;
      }

      for (int w : adj[dir][u]) {
        if(A[w] != 999) {
          temp[w] = temp[u] + 1;
          if(use == 1 || prunable(i, w, temp, A, dir)) {
            quee[1].push(w);
          } else {
            distances[dir][w][i] = temp[w];
            quee[0].push(w);
          }
          A[w] = 999;
        }
      }
    }
  }
}

bool DHighwayLabelling::prunable(int i, int u, int *temp, int *A, int dir) {

  if(landmarks.count(u) > 0) {
    highway[dir][i][landmarks[u]] = temp[u];
    return true;
  } else {
    for (int w : adj[1-dir][u]) {
      if(A[w] == 999) {
        if(temp[w] == temp[u] - 1 && distances[dir][w][i] == 999)
          return true;
      }
    }
  }
  return false;
}

int DHighwayLabelling::query(int r, int v, int dir) {

  int m = 999;
  for(int i = 0; i < K; i++) {
    m = min(m, distances_1[dir][v][i] + highway_1[dir][r][i]);
  }
  return m;
}

int DHighwayLabelling::min(int a, int b) {
  return (a < b) ? a : b;
}

void DHighwayLabelling::SelectLandmarks_HD(int topk[]) {
  std::vector<std::pair<long, int> > deg(V); long long sum = 0;
  for (int v = 0; v < V; v++) {
    deg[v] = std::make_pair(adj[0][v].size() + adj[1][v].size(), v);
    sum = sum + (adj[0][v].size() + adj[1][v].size());
  }
  std::sort(deg.rbegin(), deg.rend());
  std::cout << "Density : " << (double) (E / 2) / V << " Avg. Degree : " << (double) sum / V << " Max. Degree : " << deg[0].first << std::endl;

  for (int v = 0; v < K; v++)
    topk[v] = deg[v].second;
}

void DHighwayLabelling::RemoveLandmarks(int landmarks[]) {

  for(int dir = 0; dir < 2; dir++) {
    for(int i = 0; i < K; i++) {
      for (int v : adj[dir][landmarks[i]]) {
        adj[dir][v].erase(std::remove(adj[dir][v].begin(), adj[dir][v].end(), landmarks[i]), adj[dir][v].end());
        adj[dir][v].shrink_to_fit();
      }
      adj[dir][landmarks[i]].clear();
      adj[dir][landmarks[i]].shrink_to_fit();
    }
  }
}

int DHighwayLabelling::HC_UB_naive(int s, int t) {

  int m = 999; int i, j;
  for(i = 0; i < C[0][s]; i++) {
    for (j = 0; j < C[1][t]; j++) {
      m = min(m, min(distances[0][s][i] + highway[0][vertices[0][s][i]][vertices[1][t][j]] + distances[1][t][j], distances[0][s][i] + highway[1][vertices[0][s][i]][vertices[1][t][j]] + distances[1][t][j]));
    }
  }
  return m;
}

void DHighwayLabelling::QueryDistance(std::string pairs, std::string output) {
  std::vector<TwoLayerQueue> que; std::vector<int> dist[2];

  dist[0].resize(V, 999); dist[1].resize(V, 999);
  que.push_back(TwoLayerQueue(V)); que.push_back(TwoLayerQueue(V));

  time_querying_millisec_ = 0; int s = 0, t = 0; int total = 0;
  std::ifstream ifs(pairs); std::ofstream ofs(std::string(output) + std::to_string(K) + std::string(".txt"));
  while(ifs >> s >> t) { total++;

    double a = -GetCurrentTimeMilliSec();
    int dist_upper = HC_UB_naive(s, t);

    int res = dist_upper, dis[2] = {0, 0};
    for (int dir = 0; dir < 2; dir++){
      int v= dir == 0 ? s : t;
      que[dir].clear();
      que[dir].push(v);
      que[dir].next();
      dist[dir][v] = 0;
    }

    while (!que[0].empty() && !que[1].empty()) {
      int use = 0;
      use = (que[0].size() <= que[1].size()) ? 0 : 1;
      dis[use]++;

      if (dis[0] + dis[1] == dist_upper) {
        res = dis[0] + dis[1];
        goto LOOP_END;
      }

      while (!que[use].empty()) {
        int v = que[use].front();
        que[use].pop();

        for (int w : adj[use][v]) {

          int &src_d = dist[use][w];
          int &dst_d = dist[1 - use][w];
          if (src_d != 999) continue;
          if (dst_d != 999) {
            res = dist[use][v] + 1 + dst_d;
            goto LOOP_END;
          } else {
            que[use].push(w);
            dist[use][w] = dist[use][v] + 1;
          }
        }
      }
      que[use].next();
    }
    LOOP_END:

    a += GetCurrentTimeMilliSec();
    time_querying_millisec_ += a;

    for (int dir = 0; dir < 2; dir++) {
      for (int v : que[dir]) {
        dist[dir][v] = 999;
      }
      que[dir].clear();
    }

    ofs << s << " " << t << " " << (int) min(res, dist_upper) << "\n";
  }
  std::cout << "QueryTime[ms]: " << (double) time_querying_millisec_ / total << std::endl;

  for(int dir = 0; dir < 2; dir++) {
    for(int i = 0; i < V; i++) {
      delete [] distances[dir][i];
      delete [] vertices[dir][i];
    }
    delete [] distances[dir];
    delete [] vertices[dir];
    delete [] C[dir];

    for(int i = 0; i < K; i++)
      delete [] highway[dir][i];
    delete [] highway[dir];
  }
}

void DHighwayLabelling::loadLabelling_Full(std::string filename, int topk[]) {

  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;

  std::ifstream ifs_1(std::string(filename) + std::to_string(K) + std::string("_index"));
  std::ifstream ifs_2(std::string(filename) + std::to_string(K) + std::string("_highway"));

  for (int dir = 0; dir < 2; dir++) {

    distances[dir] = new int*[V]; distances_1[dir] = new int*[V];
    for(int i = 0; i < V; i++) {
      distances[dir][i] = new int[K]; distances_1[dir][i] = new int[K];
      for(int j = 0; j < K; j++) {
        distances[dir][i][j] = 999;
        distances_1[dir][i][j] = 999;
      }
    }

    int C, idx;
    for(int i = 0; i < V; i++) {
      ifs_1.read((char*)&C, sizeof(C));
      for(int j = 0; j < C; j++) {
        ifs_1.read((char*)&idx, sizeof(idx));
        ifs_1.read((char*)&distances[dir][i][idx], sizeof(distances[dir][i][idx]));
        distances_1[dir][i][idx] = distances[dir][i][idx];
      }
    }

    highway[dir] = new int*[K]; highway_1[dir] = new int*[K];
    for(int i = 0; i < K; i++) {
      highway[dir][i] = new int[K]; highway_1[dir][i] = new int[K];
      for(int j = 0; j < K; j++) {
        ifs_2.read((char*)&highway[dir][i][j], sizeof(highway[dir][i][j]));
        highway_1[dir][i][j] = highway[dir][i][j];
      }
    }
  }
  ifs_1.close();
  ifs_2.close();

  std::cout << LabellingSize() << std::endl;
}

void DHighwayLabelling::loadLabelling_Pruned(std::string filename) {
  std::ifstream ifs_1(std::string(filename) + std::to_string(K) + std::string("_index"));
  std::ifstream ifs_2(std::string(filename) + std::to_string(K) + std::string("_highway"));

  for (int dir = 0; dir < 2; dir++) {
    C[dir] = new int[V];
    vertices[dir] = new int*[V];
    distances[dir] = new int*[V];

    for(int i = 0; i < V; i++) {
      ifs_1.read((char*)&C[dir][i], sizeof(C[dir][i]));
      vertices[dir][i] = new int[C[dir][i]];
      distances[dir][i] = new int[C[dir][i]];
      for(int j = 0; j < C[dir][i]; j++) {
        ifs_1.read((char*)&vertices[dir][i][j], sizeof(vertices[dir][i][j]));
        ifs_1.read((char*)&distances[dir][i][j], sizeof(distances[dir][i][j]));
      }
    }

    highway[dir] = new int*[K];
    for(int i = 0; i < K; i++) {
      highway[dir][i] = new int[K];
      for(int j = 0; j < K; j++)
        ifs_2.read((char*)&highway[dir][i][j], sizeof(highway[dir][i][j]));
    }
  }
  ifs_1.close();
  ifs_2.close();
}

void DHighwayLabelling::storeLabelling(std::string filename) {
  std::ofstream ofs_1(std::string(filename) + std::to_string(K) + std::string("_index"));
  std::ofstream ofs_2(std::string(filename) + std::to_string(K) + std::string("_highway"));

  for (int dir = 0; dir < 2; dir++) {
    for(int i = 0; i < V; i++) {
      int C = 0;
      for(int j = 0; j < K; j++) {
        if(distances[dir][i][j] != 999)
          C++;
      }

      ofs_1.write((char*)&C, sizeof(C));
      for(int j = 0; j < K; j++) {
        if(distances[dir][i][j] != 999) {
          ofs_1.write((char*)&j, sizeof(j));
          ofs_1.write((char*)&distances[dir][i][j], sizeof(distances[dir][i][j]));
        }
      }
    }

    for(int i = 0; i < K; i++) {
      for(int j = 0; j < K; j++) {
        ofs_2.write((char*)&highway[dir][i][j], sizeof(highway[dir][i][j]));
      }
    }
  }
  ofs_1.close();
  ofs_2.close();
}

#endif  // DHIGHWAY_LABELING_H_
