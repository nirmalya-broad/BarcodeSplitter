#ifndef BKNODE_H
#define BKNODE_H
/*
    BKNode.h

    By Stephen Holiday 2011
    http://stephenholiday.com
    (Exception, Distance Algorithm by Anders Sewerin Johansen)

    The code is under the [Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0) license.

    This a C++ implementation of a BKTree (Burkhard-Keller).
    Essentially it allows searching for an object in a metric space.

    One common use is for fuzzy string matching, like in a spell checker.
    The search is performed by looking at the distance of the test with the current
    node and moving according to the distance of the test to the children.

    This technique is faster than brute force as it does not need to look at every
    possible node in the space.

    Implemented according to this post:
    http://blog.notdot.net/2007/4/Damn-Cool-Algorithms-Part-1-BK-Trees

*/

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>

// Functor for deleting pointers in map.
template<class A, class B>
struct DeleteMapFntor
{
    // Overloaded () operator.
    // This will be called by for_each() function.
    bool operator()(std::pair<A,B> x) const
    {
        // Assuming the second item of map is to be
        // deleted. Change as you wish.
        delete x.second;
        return true;
    }
};

template <typename T>
class BKNode {
    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
        
    protected:
        T value;
        std::map<int,BKNode<T>*> children;
    public:
        BKNode(const T &);
        BKNode();
        BKNode(T, std::map<int,BKNode<T> >);
        ~BKNode();
        
        // Accessors
        T get() const;

        int distance(const BKNode<std::string> &) const;
        std::vector<T> find(const T &, const int);
        
        // Mutators
        void insert(const T &);
        
    template <typename TT>
    friend std::ostream &operator << ( std::ostream &, const BKNode<TT> & );
    
};

template <typename T>
template<class Archive>
void BKNode<T>::serialize(Archive & ar, const unsigned int version)
{
    ar & value;
    ar & children;
}

template <typename T>
BKNode<T>::BKNode() {
}

template <typename T>
BKNode<T>::BKNode(const T &newValue) {
    value=newValue;
    children=std::map<int,BKNode<T>*>();
}

template <typename T>
BKNode<T>::BKNode(T v, std::map<int,BKNode<T> > c) {
    value=v;
    children=c;
}

template <typename T>
BKNode<T>::~BKNode() {
    for_each( children.begin(),children.end(),DeleteMapFntor<int,BKNode<T>*>());
}

template <typename T>
T BKNode<T>::get() const{
    return value;
}

template <typename T>
void BKNode<T>::insert(const T &rhs) {
    int dist = distance(rhs);
    
    // If the distance is zero, they're the same objects so don't insert
    if (dist!=0) {
        
        // If there is already a child with the same distance, insert it on that child
        if (children.count(dist)) {
            children[dist]->insert(rhs);
        } else {
            BKNode<T> *new_node = new BKNode<T>(rhs);
            children[dist]=new_node;
        }
    }
}


template <typename T>
std::vector<T> BKNode<T>::find(const T &rhs, const int threshold) {
    int dist = distance(rhs);
    
    std::vector<T> results=std::vector<T>();
    // If the current node is in the threshold, return it
    if (dist<=threshold)
        results.push_back(value);
    
    int dmin=dist-threshold;
    int dmax=dist+threshold;
    
    for (int i=dmin; i<=dmax; i++) {
        if (children.count(i)) {
            std::vector<T> partial= std::vector<T>(children[i]->find(rhs,threshold));
            
            results.insert(results.end(), partial.begin(), partial.end());
        }
    }
    
    return results;
}


// By Nirmalya: Since for allseq we only need substitution, we shall
// consider L1 norm. We shall not consider insertion, deletion and
// transposition.

template <>
int BKNode<std::string>::distance(const BKNode<std::string>& rhs) const{
    std::string source=value;
    std::string target=rhs.value;
  
    /*
        From http://www.merriampark.com/ldcpp.htm
        By Anders Sewerin Johansen
    */
      // Step 1

    const int n = source.length();
    const int m = target.length();
    if (n == 0) {
        throw std::invalid_argument("Length of source is zero.");
    }
    if (m == 0) {
        throw std::invalid_argument("Length of target is zero.");
    }

	if (m !=n ){ 
		std::string msg = "Source and target have different length.\n" 
			"The size of the barcode from file1 does not match with\n" 
			" one from the dictionary.\n" + source + "\n" + target + "\n";
		throw std::invalid_argument(msg);
	}
	int ldist = 0;

	for (int j = 0; j < n; j++) {
		if (source[j] != target[j])	{
			ldist++;
		}
	}

	return ldist;
}

template <typename TT>
std::ostream &operator << ( std::ostream &out, const BKNode<TT> &node ) {
    out << "BKNode<this: " << node.value << " >";
    
    return out;
}

#endif
