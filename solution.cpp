#pragma GCC optimize("Ofast", "unroll-loops", "omit-frame-pointer", "inline")
#pragma GCC option("arch=native", "tune=native", "no-zero-upper")
#pragma GCC target("rdrnd", "popcnt", "avx", "bmi2")

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include <utility>
#include <sys/time.h>

using namespace std;

/**
 * Auto-generated code below aims at helping you parse VERSION ============
 * the standard input according to the problem statement.
 **/

#define CRYSTAL 2
#define EGGES 1
#define ALL 10
int STRENGTH;
typedef struct t_cell {
  int type;              // 0 for empty, 1 for eggs, 2 for crystal
  int initial_resources; // the initial amount of eggs/crystals on this cell

  int neigh[6];

  int resources; // the current amount of eggs/crystals on this cell
  int opp_ants;  // the amount of opponent ants on this cell

  
  int my_ants;   // the amount of your ants on this cell
  int beacons = 0;
  int indx;
  int wiggleRoom;

} t_cell;

typedef struct t_single_data {
  
  int all_my_ants;
  int my_intial_ants;
  int all_opp_ants;
  int all_intial_crystal;
  int all_intial_egges;
  int turn_crystal;
  int number_of_cells;
  int number_of_bases;
  int myScore;
  int oppScore;
  int topScor;
  int **hash_tab;

} t_single_data;
/******************************************/

struct AntBeaconPair {
  t_cell *ant;
  t_cell *beacon;
  // ... other members
};

struct AntAllocation {
  int antIndex;
  int beaconIndex;
  int count;
  // ... other members
};

int getDistance(vector<t_cell> vec, t_single_data data, int start, int target) {
  vector<bool> visited(data.number_of_cells, false);
  vector<int> parent(data.number_of_cells, -1);
  std::queue<int> q;
  q.push(start);
  visited[start] = true;
  bool pathFound = false;
  while (!q.empty()) {
    int current = q.front();
    q.pop();
    if (current == target) {
      pathFound = true;
      break;
    }
    for (auto index : vec[current].neigh) {
      if (index != -1 && !visited[index]) {
        q.push(index);
        visited[index] = true;
        parent[index] = current;
      }
    }
  }
  vector<int> path;
  if (pathFound) {
    int current = target;
    while (current != start) {
      path.push_back(current);
      current = parent[current];
    }
    path.push_back(start);
  }
  if (path.size() != 0)
    return (path.size() - 1);
  return (-1);
}
/************************************************************************************************************/
std::vector<AntAllocation> innerAllocateAnts(vector<t_cell> vec, t_single_data data, std::vector<t_cell> antCells,std::vector<t_cell> beaconCells)
{
  // cerr<<"innerAllocateAnts_1"<<endl;
  std::vector<AntAllocation> allocations;
  long long antSum = 0;
  for (const t_cell &cell : antCells) {
    antSum += cell.my_ants;
  }

  long long beaconSum = 0;
  for (const t_cell &cell : beaconCells) {
    beaconSum += cell.beacons;
  }

  double scalingFactor = (double)(antSum)/beaconSum;

    for(int i = 0; i < beaconCells.size(); i++)
    {
        long long highBeaconValue = static_cast<long long>(ceil(beaconCells[i].beacons * scalingFactor));
        long long lowBeaconValue = static_cast<long long>(beaconCells[i].beacons * scalingFactor);
        beaconCells[i].beacons = std::max(1LL, lowBeaconValue);
        beaconCells[i].wiggleRoom = highBeaconValue - beaconCells[i].beacons;
    }




  std::vector<AntBeaconPair> allPairs;

  for (t_cell &antCell : antCells) {
    for (t_cell &beaconCell : beaconCells) {
      AntBeaconPair pair{&antCell, &beaconCell};
      allPairs.push_back(pair);
    }
  }
 
  std::sort(allPairs.begin(), allPairs.end(),
            [data](AntBeaconPair &pair1, AntBeaconPair &pair2) {
              double distance1 = data.hash_tab[pair1.ant->indx][pair1.beacon->indx];
              double distance2 = data.hash_tab[pair2.ant->indx][pair2.beacon->indx];
              if (distance1 != distance2)
                return distance1 < distance2;
              if (pair1.ant->indx != pair2.ant->indx)
                return pair1.ant->indx < pair2.ant->indx;
              return pair1.beacon->indx < pair2.beacon->indx;
            });

  bool stragglers = false;
// cerr<<"innerAllocateAnts_2"<<endl;
  int i = 0;
  while (antSum > 0) {
    for (AntBeaconPair &pair : allPairs) {
      t_cell *antCell = pair.ant;
      long long antCount = antCell->my_ants;
      t_cell *beaconCell = pair.beacon;
      long long beaconCount = beaconCell->beacons;
      long long wiggleRoom = beaconCell->wiggleRoom;

      int maxAlloc = (stragglers ? std::min(antCount, beaconCount + wiggleRoom): std::min((antCount), (beaconCount)));

      if (maxAlloc > 0) {
        allocations.push_back(AntAllocation{antCell->indx, beaconCell->indx, maxAlloc});
        

        antSum -= maxAlloc;
        antCell->my_ants -= maxAlloc;
        if (!stragglers) {
          beaconCell->beacons -= maxAlloc;
        } else {
          beaconCell->beacons -= (maxAlloc - wiggleRoom);
          beaconCell->wiggleRoom = 0;
        }
      }
    }
    stragglers = true;
  }
  // cerr<<"innerAllocateAnts_3"<<endl;
  return allocations;
}
/************************************************************************************************************/
int  ft_get_next_move(vector<t_cell> &vec, t_single_data data, int a, int b)
{
    int dist = data.hash_tab[a][b];
    vector<int> neibr_celles;
    for(auto neib : vec[a].neigh)
    {
        if (neib != -1 && dist > data.hash_tab[neib][b])
           neibr_celles.push_back(neib);
    }
    if (neibr_celles.empty())
        return (a);
    if (neibr_celles.size() == 1)
        return (neibr_celles[0]);
    
    if(vec[neibr_celles[0]].my_ants != vec[neibr_celles[1]].my_ants)
        return ((vec[neibr_celles[0]].my_ants > vec[neibr_celles[1]].my_ants) ? neibr_celles[0] : neibr_celles[1]);
    
    if (vec[neibr_celles[0]].beacons != vec[neibr_celles[1]].beacons)
        return ((vec[neibr_celles[0]].beacons > vec[neibr_celles[1]].beacons) ? neibr_celles[0] : neibr_celles[1]);
    
    return (min(neibr_celles[0], neibr_celles[1]));
   
}
/*******************************************************/
std::vector<t_cell> ft_get_antCells(vector<t_cell> vec) {
  vector<t_cell> antCells;
  for (auto cell : vec)
    if (cell.my_ants)
      antCells.push_back(cell);

  return (antCells);
}
/*********************************************************/
std::vector<t_cell> ft_get_beaconCells(vector<t_cell> vec) {
  vector<t_cell> beaconCells;
  for (auto cell : vec)
    if (cell.beacons)
      beaconCells.push_back(cell);

  return (beaconCells);
}
/*********************************************************/
int ft_do_beacons(vector<t_cell> vec) {
  int ret = 0;

  for (auto indx : vec)
    if (indx.beacons)
      cout << "BEACON " << indx.indx << " " << indx.beacons << ";";

  return (ret);
}
/*********************************************************************************************************************/
/*********************************************************************************************************************/
int ft_get_min(vector<t_cell> vec, vector<int> path, int red_blu)
{
    int min = 9999;
    for(auto index : path)
    {
        if ((vec[index].my_ants*(red_blu == 1) + vec[index].opp_ants*(red_blu == 0)) < min)
            min = (vec[index].my_ants*(red_blu == 1) + vec[index].opp_ants*(red_blu == 0));
    }
    if (min == 9999)
        min = 0;
    return(min);
}
/*****************************************************************************************/
vector<int> ft_stronget_path(vector<t_cell> vec, t_single_data data, int cell, vector<int> bases, int red_blu)
{
    vector<bool> visited(data.number_of_cells, false);
    vector<int> parent(data.number_of_cells, -1);
    std::queue<int> q;
    q.push(cell);
    visited[cell] = true;
    bool pathFound = false;
    int base;
    while(!q.empty())
    {
        int current = q.front();
        q.pop();
        if (current == bases[0] || current == bases[1]) {
            base = current == bases[0] ? bases[0] : bases[1];
            pathFound = true;
            break;
        }
        for(auto index : vec[current].neigh)
        {
            if(index != -1 && !visited[index] && (vec[index].my_ants*(red_blu == 1) || vec[index].opp_ants*(red_blu == 0)))
            {
                q.push(index);
                visited[index] = true;
                parent[index] = current;
            }
        }

    }
    vector<int> path;
    if (pathFound)
    {
        int current = base;
        while(current != cell)
        {
            path.push_back(current);
            current = parent[current];
        }
        path.push_back(cell);
    }
    return(path);
}

/*******/
std::pair<int, int> ft_get_scor(vector<t_cell> vec, t_single_data data, vector<int> my_base, vector<int> opp_base)
{
    vector<int> chien;
    int min = 987654321;
    int win_egges = 0;
    int win_crystal = 0;
    std::pair<int, int> myPair;
    for (auto cell : vec)
    {
        if (cell.my_ants && cell.resources && cell.type == 1)
        {
            chien = ft_stronget_path(vec, data, cell.indx, my_base, 1);
            for (auto chien_indx : chien)
            {
                if (vec[chien_indx].my_ants < min)
                {
                    min = vec[chien_indx].my_ants;
                }
                if (vec[chien_indx].opp_ants)
                {
                    if (min < ft_get_min(vec, ft_stronget_path(vec, data, cell.indx, opp_base, 0), 0))
                        min = 0;
                }
                if (min == 0)
                    break;
            }
            if (chien.size() && min)
            {
            win_egges += (min <= cell.resources ? min : cell.resources);
            }
        }
    }
    vec[my_base[0]].my_ants += win_egges;
    if (my_base[1] != -10)
        vec[my_base[1]].my_ants += win_egges;
    min = 987654321;
    for (auto cell : vec)
    {
        if (cell.my_ants && cell.resources && cell.type == 2)
        {
            chien = ft_stronget_path(vec, data, cell.indx, my_base, 1);
            for (auto chien_indx : chien)
            {
                if (vec[chien_indx].my_ants < min)
                {
                    min = vec[chien_indx].my_ants;
                }
                if (vec[chien_indx].opp_ants)
                {
                    if (min < ft_get_min(vec, ft_stronget_path(vec, data, cell.indx, opp_base, 0), 0))
                        min = 0;
                }
                if (min == 0)
                    break;
            }
            if (chien.size() && min)
            {
            win_crystal += (min <= cell.resources ? min : cell.resources);
            }
        }
    }
    myPair.first = win_egges;
    myPair.second = win_crystal;

    return (myPair);
}
/*********************************************************************************************************************/
/*********************************************************************************************************************/
vector<t_cell> ft_change_state(vector<t_cell> vec, t_single_data data, std::vector<AntAllocation> ant_alocation)
{
    vector<t_cell> vec_copy1 = vec;
    vector<t_cell> vec_copy2 = vec;
    for(int i = 0; i < ant_alocation.size() ; i++)
    {
        int next_indx = ft_get_next_move(vec_copy1, data, ant_alocation[i].antIndex, ant_alocation[i].beaconIndex);
        vec_copy2[ant_alocation[i].antIndex].my_ants -= ant_alocation[i].count;
        vec_copy2[next_indx].my_ants += ant_alocation[i].count;
    }
    return (vec_copy2);
} 
/*********************************************************************************************************************/
/*********                                              ft_random                                               ******/
/*********************************************************************************************************************/
int ft_lala(vector<t_cell> vec, int indx)
{
    if (vec[indx].my_ants)
        return (1);
    for (auto neib : vec[indx].neigh)
        if (vec[neib].my_ants)
            return (1);
    return(0);
}
void ft_random(vector<t_cell> &vec)
{
    int beaeconed_number = 0;
    for (int i = 0; i < vec.size() ; i++)
    {
        if (vec[i].beacons)
        {
          int a = rand() % 80 + 10;
            vec[i].beacons = a;
            beaeconed_number++;
        }
    }
    
}
/********/
void ft_random_beacon_set(vector<t_cell> &vec, int number_of_cells)
{
    int num_sets = 0;
    for (int i = 0; i < vec.size() ; i++)
    {
        if (num_sets < 10)
        {
          int indx = rand() % number_of_cells;
            vec[indx].beacons = rand() % 80 + 10;;
            num_sets++;
        }
    }
    
}
/*********************************************************************************************************************/
/*********************************************************************************************************************/
/******************************************-BFS-***********************************************/
vector<int> ft_BFS_closetResource(vector<t_cell> vec, t_single_data data, int start, int order, int type = 2)
{
    vector<bool> visited(data.number_of_cells, false);
    vector<int> parent(data.number_of_cells, -1);
    std::queue<int> q;
    q.push(start);
    visited[start] = true;
    bool pathFound = false;
    int next_target = 1;
    int target;
    while(!q.empty())
    {
        int current = q.front();
        q.pop();
        if (vec[current].resources && (vec[current].type == type || type == 10) && !vec[current].beacons) {
            if (next_target == order)
            {
                pathFound = true;
                target = current;
                break;
            }
            next_target++;
        }
        for(auto index : vec[current].neigh)
        {
            if(index != -1 && !visited[index])
            {
                q.push(index);
                visited[index] = true;
                parent[index] = current;
            }
        }

    }
    vector<int> path;
    if (pathFound)
    {
        int current = target;
        while(current != start)
        {
            path.push_back(current);
            current = parent[current];
        }
        path.push_back(start);
    }
    return(path);
}
/*************************************************************************************************/
/*************************************************************************************************/

int ft_how_beacon_around(vector<t_cell> &vec, int cell)
{
    int ret = 0;
    for(auto indx : vec[cell].neigh)
        if (vec[indx].beacons && indx != -1)
            ret++;
    for(int i = 0; i < 5; i++)
            if (vec[vec[cell].neigh[i]].beacons && vec[vec[cell].neigh[i+1]].beacons && vec[cell].neigh[i] != -1 && vec[cell].neigh[i+1] != -1)
            ret -= 1;
    if (vec[vec[cell].neigh[0]].beacons && vec[vec[cell].neigh[5]].beacons && vec[cell].neigh[0] != -1 && vec[cell].neigh[5] != -1)
            ret -= 1;
    return (ret);
}


int ft_hallf_egges(vector<t_cell> vec, t_single_data data)
{
  int egges_touched = 0;
  for (auto cell: vec)
    if(cell.beacons && cell.type == 1)
      egges_touched += cell.resources;
  egges_touched += (data.all_my_ants - data.my_intial_ants)/data.number_of_bases;

  return (egges_touched);
}

int ft_hallf_crystals(vector<t_cell> vec, t_single_data data)
{
  int crystals_touched = 0;
  for (auto cell: vec)
    if(cell.beacons && cell.type == 2)
      crystals_touched += cell.resources;
  crystals_touched += data.myScore;

  return (crystals_touched);
}
/*************************************************************************************************/
/*************************************************************************************************/
/*************************************************************************************************/
/*************************************************************************************************/

void ft_full_beacons(vector<t_cell> &vec, t_single_data data, vector<int> my_bases, int STRENGTH)
{
    int all_beacons_cells = 0;
    int all_my_ants = data.all_my_ants;
    for (int i = 0; i < vec.size(); i++)
    {
        int tru = 0;
        for (auto cell: vec)
        {
            if (ft_how_beacon_around(vec, cell.indx) < 2 && vec[cell.indx].beacons && cell.indx != my_bases[0] && cell.indx != my_bases[1] && !cell.resources)
            {
                vec[cell.indx].beacons = 0;
                tru = 1;
            }
        }
        if (tru == 0)
            break;
    }
    for(auto cell :vec)
    {
        if (vec[cell.indx].beacons)
        {
            all_beacons_cells++;
        }
    }
    
    vector<int> short_path;
    while (all_my_ants/(STRENGTH) >= all_beacons_cells && ft_hallf_egges(vec, data) < 2*data.all_intial_egges/3 && data.oppScore < data.topScor/3)/////////////////egges loop
    {
        int short_size = 2500;
        vector<int> tmp;
        short_path.clear();
        for(auto cell :vec)
        {
            if (cell.beacons)
            {
                tmp = ft_BFS_closetResource( vec, data, cell.indx, 1, EGGES);
                if (tmp.size() < short_size && tmp.size())
                {
                    short_path = tmp;
                    short_size = tmp.size();
                }
            }
            
        }
        if(short_path.size())
            short_path.pop_back();
        if (all_my_ants/(STRENGTH) < (all_beacons_cells + short_path.size()) || short_path.empty())
        {
            break;
        }
        all_beacons_cells += short_path.size();
        for(auto indx : short_path)
        {
            vec[indx].beacons = STRENGTH;
        }
    }

    while (all_my_ants/(STRENGTH) >= all_beacons_cells && ft_hallf_crystals(vec, data) < data.topScor || all_beacons_cells == data.number_of_bases && data.oppScore >= data.topScor/4)/////////////////crystal loop
    {
        int short_size = 2500;
        vector<int> tmp;
        short_path.clear();
        for(auto cell :vec)
        {
            if (cell.beacons)
            {
                tmp = ft_BFS_closetResource( vec, data, cell.indx, 1, CRYSTAL);
                if (tmp.size() < short_size && tmp.size())
                {
                    short_path = tmp;
                    short_size = tmp.size();
                }
            }
            
        }
        if(short_path.size())
            short_path.pop_back();
        if (all_my_ants/(STRENGTH) < (all_beacons_cells + short_path.size()) || short_path.empty())
        {
            break;
        }
        all_beacons_cells += short_path.size();
        for(auto indx : short_path)
        {
            vec[indx].beacons = STRENGTH;
        }
    }
}
/*********************************************************************************************************************/
/*********************************************************************************************************************/
int main() {
  srand(time(NULL));
  int number_of_cells; // amount of hexagonal cells in this map
  cin >> number_of_cells;
  cin.ignore();
  vector<t_cell> vec(number_of_cells);
  t_single_data data;
  data.number_of_cells = number_of_cells;
  for (int i = 0; i < number_of_cells; i++) {
    t_cell cell;
    int type;              // 0 for empty, 1 for eggs, 2 for crystal
    int initial_resources; // the initial amount of eggs/crystals on this cell
    int neigh_0; // the index of the neighbouring cell for each direction
    int neigh_1;
    int neigh_2;
    int neigh_3;
    int neigh_4;
    int neigh_5;
    cin >> type >> initial_resources >> neigh_0 >> neigh_1 >> neigh_2 >>neigh_3 >> neigh_4 >> neigh_5;cin.ignore();
    vec[i].type = type; // 0 for empty, 1 for eggs, 2 for crystal
    vec[i].initial_resources = initial_resources; // the initial amount of eggs/crystals on this cell
    vec[i].neigh[0] = neigh_0; // the index of the neighbouring cell for each direction
    vec[i].neigh[1] = neigh_1;
    vec[i].neigh[2] = neigh_2;
    vec[i].neigh[3] = neigh_3;
    vec[i].neigh[4] = neigh_4;
    vec[i].neigh[5] = neigh_5;
    vec[i].indx = i;
  }
  int number_of_bases;
  cin >> number_of_bases;cin.ignore();
  vector<int> my_bases(2, -10);
  vector<int> opp_bases(2, -10);
  for (int i = 0; i < number_of_bases; i++) {
    int my_base_index;
    cin >> my_base_index;cin.ignore();
    my_bases[i] = my_base_index;
  }
  for (int i = 0; i < number_of_bases; i++) {
    int opp_base_index;
    cin >> opp_base_index;cin.ignore();
    opp_bases[i] = opp_base_index;
  }

  // game loop
  int my_base_indx;
  int opp_base_indx;
  int get_base = 0;
  int turn = 0;
  data.all_intial_crystal = 0;
  data.all_intial_egges = 0;
  data.my_intial_ants = 0;
  data.number_of_bases = number_of_bases;

int **hash_tab = (int **)malloc(sizeof(int*)*number_of_cells);
for(int i = 0; i < number_of_cells ;i++)
{
    hash_tab[i] = (int *)malloc(sizeof(int)*number_of_cells);
    for(int j = 0; j < number_of_cells; j++)
    {
        hash_tab[i][j] = getDistance(vec, data, i, j);
    }
}
data.hash_tab = hash_tab;



  while (1) {
    turn++;
    data.all_my_ants = 0;
    data.all_opp_ants = 0;
    data.turn_crystal = 0;
    int action = 0;
    cin >> data.myScore >> data.oppScore;cin.ignore();
    for (int i = 0; i < number_of_cells; i++) {
      int resources; // the current amount of eggs/crystals on this cell
      int my_ants;   // the amount of your ants on this cell
      int opp_ants;  // the amount of opponent ants on this cell
      cin >> resources >> my_ants >> opp_ants;cin.ignore();

      if (turn == 1)
      {
        if (vec[i].type == 2)
          data.all_intial_crystal += resources;
        if (vec[i].type == 1)
          data.all_intial_egges += resources;
        data.my_intial_ants += my_ants;
      }

      vec[i].beacons = 0;
      vec[i].resources = resources;
      vec[i].type *= (resources != 0);
      vec[i].my_ants = my_ants;
      vec[i].opp_ants = opp_ants;
      data.all_my_ants += my_ants;
      data.all_opp_ants += opp_ants;
      if (vec[i].type == 2)
        data.turn_crystal += resources;
    }
    data.topScor = data.all_intial_crystal/2;
    struct timeval start, end;
    gettimeofday(&start, nullptr);
    long long startMs = start.tv_sec * 1000 + start.tv_usec / 1000;
    // vec[0].all_my_ants /= number_of_bases;
    /****************************/

    STRENGTH = 3;
    if(ft_BFS_closetResource( vec, data, my_bases[0], 1, EGGES).size() > data.all_my_ants/STRENGTH)
        STRENGTH = 2;
    if(ft_BFS_closetResource( vec, data, my_bases[0], 1, EGGES).size() > data.all_my_ants/2)
        STRENGTH = 1; 
    vec[my_bases[0]].beacons = STRENGTH;
    if (number_of_bases == 2)
        vec[my_bases[1]].beacons = STRENGTH;
    /*******************************************************************************/
    int i = 0;
        ft_full_beacons(vec, data, my_bases, STRENGTH);
        std::vector<AntAllocation> ant_alocation = innerAllocateAnts(vec, data, ft_get_antCells(vec), ft_get_beaconCells(vec));
        vector<t_cell> vec_copy = vec;
        std::pair<int, int> my_scor = ft_get_scor(ft_change_state(vec, data, ant_alocation), data, my_bases, opp_bases);
        while(++i)
        {

              ft_random(vec);

              std::vector<AntAllocation> ant_alocation = innerAllocateAnts(vec, data, ft_get_antCells(vec), ft_get_beaconCells(vec));

              std::pair<int, int> tmp_scor = ft_get_scor(ft_change_state(vec, data, ant_alocation), data, my_bases, opp_bases);


              if (tmp_scor.first*2 + tmp_scor.second > my_scor.first*2 + my_scor.second)
              {
                  vec_copy = vec;
                  my_scor = tmp_scor;
              }
              
            gettimeofday(&end, nullptr);
            long long endMs = end.tv_sec * 1000 + end.tv_usec / 1000;
            if(endMs - startMs >= 98 && turn != 1)
            {
               
                break;
            }
            if(endMs - startMs >= 990 && turn == 1)
            {
                break;
            }
          
             
        }
        //some code missed here
    cout << "WAIT;" << endl;
  }
}