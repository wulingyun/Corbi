#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "Corbi.h"

#include "FibHeap.h"
#include <stdlib.h>
//#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <cstdio>
using   namespace   std;

class HeapNode: public FibHeapNode
{
	double m_label;
	int m_index;

public:
	HeapNode() : FibHeapNode() { m_label = 0; };

	virtual void operator =(FibHeapNode& RHS);
	virtual int  operator ==(FibHeapNode& RHS);
	virtual int  operator <(FibHeapNode& RHS);

	virtual void operator =(double key);

	double getKeyValue() { return m_label; }
	void setKeyValue(double key) { m_label = key; }

	int getIndexValue() { return m_index; }
	void setIndexValue(int i) { m_index = i; }
};

void HeapNode::operator =(double label)
{
	HeapNode temp;
	temp.m_label = m_label = label;
	FHN_assign(temp);
}

void HeapNode::operator =(FibHeapNode& RHS)
{
	FHN_assign(RHS);
	m_label = ((HeapNode&) RHS).m_label;
}

int  HeapNode::operator ==(FibHeapNode& RHS)
{
	if (FHN_compare(RHS)) return 0;
	return m_label == ((HeapNode&) RHS).m_label ? 1 : 0;
}

int  HeapNode::operator <(FibHeapNode& RHS)
{
	int X;

	if ((X=FHN_compare(RHS)) != 0) return X < 0 ? 1 : 0;
	return m_label < ((HeapNode&) RHS).m_label ? 1 : 0;
}


void onenoded(double *W,int d1,int s,double *D){
	// ** the pointer of the weight matrix
	// d1 the dimension of the matrix
	// s the source node 
	// * the one row pointer of the distance matrix 

    HeapNode *label, *min_label, temp_label;
	FibHeap heap;
	double  current_value,d;	
	int *backtrack;
	int i, j;
	bool *label_free;

	label_free = new bool [d1];
	backtrack = new int [d1];
	label = new HeapNode [d1];


	// Dijkstra's algorithm
	// initial the nodes
		for (i=0; i<d1; i++) {
			if(i!=s){
		    label[i].setKeyValue(HUGE_VAL);
			label[i].setIndexValue(i);
			heap.insert(&label[i]);
			label_free[i] = true;
			backtrack[i] = -1;}
			else{
		    label[i].setKeyValue(0);
			label[i].setIndexValue(i);
			heap.insert(&label[i]);
			label_free[i] = true;
			backtrack[i] = -1;			
			}
		}

    // do while
		while(true){	
			min_label = (HeapNode *)heap.extractMin();				
			if (min_label == NULL || min_label->getKeyValue() >= HUGE_VAL) break;
			i = min_label->getIndexValue();
			current_value = min_label->getKeyValue();
			label_free[i] = false;
			D[i+s*d1]=min_label->getKeyValue();
			
            for (j=0;j<d1;j++){
				if(W[i*d1+j]!=0 && label_free[j]){ // i ,j is neighbor
					d = current_value + W[i*d1+j];
					if (d < label[j].getKeyValue()) {
							temp_label = label[j];
							temp_label.setKeyValue(d);
							heap.decreaseKey(&label[j], temp_label);
							backtrack[j] = i;
						}
				}
			}
		}

	delete[] label;
	delete[] backtrack;
	delete[] label_free;
}

SEXP dis(SEXP _W, SEXP _d1)
{
	PROTECT(_W = AS_NUMERIC(_W));
	double *W = NUMERIC_POINTER(_W);
	int d1 = INTEGER_POINTER(AS_INTEGER(_d1))[0];

	SEXP _D;
	PROTECT(_D = NEW_NUMERIC(d1*d1));
	setDim2(_D, d1, d1);
	double *D = NUMERIC_POINTER(_D);

	int i,j;

	
	for (i=0;i<d1;i++){for(j=0;j<d1;j++){D[i+j*d1]=0;}}
	
	for (i=0;i<d1;i++){
	onenoded(W,d1,i,D);
	}

	UNPROTECT(2);
	return (_D);
}