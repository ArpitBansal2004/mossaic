/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <deque>
#include <iostream>

using namespace std;


//helper functions

//copy function for copy constructor and assignment operator
template<int Dim>
void KDTree<Dim>::copy(KDTreeNode *lhs, KDTreeNode *rhs) 
{
  lhs = new KDTreeNode();
  lhs -> point = rhs -> point;
  copy(lhs -> left, rhs -> left);
  copy(lhs -> right, rhs -> right);
}

//clear function to help in destructor and assignment operator
template<int Dim>
void KDTree<Dim>::clear(KDTreeNode *root) 
{
  //check null
  if (root == NULL) {
    return;
  }
  //call on left
  if (root -> left != NULL) {
    clear(root->left);
  }

  //call on right
  if (root -> right != NULL) {
    clear(root -> right);
  }

  //handle passed node
  delete root;  
  root = NULL;
}


template <int Dim>
void KDTree<Dim>::swap(Point<Dim> & point_1, Point<Dim> & point_2) 
{
  Point<Dim> temp = point_1;
  point_1 = point_2;
  point_2 = temp;
}

template <typename RandIter, typename Comparator>
RandIter partition(RandIter start, RandIter end, RandIter pivotIndex, Comparator cmp) 
{
  auto pivot = *pivotIndex;

  //Move pivot to the end
  std::swap(*pivotIndex, *end);
  auto store = start;
  while (start != end) {
    if (cmp(*start, pivot)) {
      std::swap(*store, *start);
      store++;
    }
    ++start;
  }
  // Put pivot element in its final position
  std::swap(*end, *store);
  return store;
}

template <int Dim>
void KDTree<Dim>::buildTree( vector<Point<Dim>>& newPoints, KDTreeNode*& curr, int dim) 
{
    //base case
    if (newPoints.empty()) {
        return;
    }

    auto cmp = [dim](const auto& arg1, const auto& arg2){
      return smallerDimVal(arg1, arg2, dim);
    };

    //calculate median index
    int median_index = floor((newPoints.size() - 1)/2);

    //Find the median element 
    nth_element(newPoints.begin(), newPoints.begin() + median_index, newPoints.end(), cmp);
    
    curr = new KDTreeNode(newPoints[median_index]);

    //create vectors for left and right subtrees
    vector<Point<Dim>> left(newPoints.begin(), newPoints.begin() + median_index);
    vector<Point<Dim>> right(newPoints.begin() + median_index + 1, newPoints.end());

    //use recursion
    buildTree(left, curr -> left, (dim + 1) % Dim);
    buildTree(right, curr -> right, (dim + 1) % Dim);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query, int dim, KDTreeNode* node) const
{
  //Handle lead nodes
  if (node -> left == NULL && node -> right == NULL) {
    return node->point;
  }
  
  //initialize curr_best and closest
  Point<Dim> curr_best = node -> point;
  Point<Dim> closest = curr_best;     
  

  //Decide whether to explore the right or left
  bool left_side = smallerDimVal(query, curr_best, dim); 

  // Recursively search the appropriate subtree
  if (node -> left != NULL && left_side) {
    closest = findNearestNeighbor(query, (dim + 1) % Dim, node -> left);
  } else if (node -> right != NULL && (! left_side)) {
    closest = findNearestNeighbor(query, (dim + 1) % Dim, node -> right);
  }

  // Update curr_best if closer than closest
  if (shouldReplace(query, curr_best, closest)) {
    curr_best = closest;
  }

  // Calculate the square of the distance along the current dimension
  double dimension_radius = (query[dim] - node -> point[dim]) * (query[dim] - node -> point[dim]);

  // Calculate the distance between query and closest
  double distance = 0;
  for (unsigned int i = 0; i < Dim; ++i) {
    distance += (query[i] - closest[i]) * (query[i] - closest[i]);
  }
  
  //Explore other subtrees along the dimension if dimension_radius is less than equal to distance
  if (dimension_radius <= distance) {
    if (node -> left != NULL && left_side == false) {
      closest = findNearestNeighbor(query, (dim + 1) % Dim, node -> left);
      if (shouldReplace(query, curr_best, closest)) {
        curr_best = closest;
      }
    } else if (node -> right != NULL && left_side == true) {
      closest = findNearestNeighbor(query, (dim + 1) % Dim, node -> right);
      if (shouldReplace(query, curr_best, closest)) {
        curr_best = closest;
      }
    }
  }

  //return
  return curr_best; // Return the nearest neighbor found so far
}


template <int Dim>
bool smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim)
{
  //check for correct dimension
  if (curDim < 0 || curDim > Dim) {
    std::cout << "incorrect dimension passed into smallerDimVal" << std::endl;
  }

  //address the tie case
  if (first[curDim] == second[curDim]) {
    return first < second;
  }

  //generic
  return first[curDim] < second[curDim];
}

template <int Dim>
bool shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential)
{

    //calculate the distance from  currentBest
    unsigned int dist_currentBest = 0;
    for (int dim_tr = 0; dim_tr < Dim; ++dim_tr) {
      dist_currentBest += ((currentBest[dim_tr] - target[dim_tr]) * (currentBest[dim_tr] - target[dim_tr]));
    }

    //caluclate the distance from potential
    unsigned int dist_potential = 0;
    for (int dim_tr = 0; dim_tr < Dim; ++dim_tr) {
      dist_potential += ((potential[dim_tr] - target[dim_tr]) * (potential[dim_tr] - target[dim_tr]));
    }

    //address the tie case (distances are equal)
    if (dist_currentBest == dist_potential) {
      return potential < currentBest;
    }

    //generic case
    return dist_potential < dist_currentBest;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    vector<Point<Dim>> newPoint_copy = newPoints;
    buildTree(newPoint_copy, root, 0);
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) 
{
  this -> size = other -> size;
  copy (this -> root, other -> root);
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) 
{
  clear(root);
  copy(root, rhs);
  this -> size = rhs -> size;
  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() 
{
  clear(root);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    return findNearestNeighbor(query, 0, root);
}

template <typename RandIter, typename Comparator>
void select(RandIter start, RandIter end, RandIter k, Comparator cmp)
{   
    //address start == end case
    if (start == end) {
      return;
    }
    auto in = start;
    auto fin = end - 1;

    //sort first
    auto pivotIndex = partition(in, fin, k, cmp);

    if (k == pivotIndex) {
      return;
    } else if (k < pivotIndex) {
      select(start, pivotIndex, k, cmp);
    } else {
      select(pivotIndex + 1, end, k, cmp);
    }
}