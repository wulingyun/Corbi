#include "Corbi.h"
#include "FibHeap.h"

class HeapNode: public FibHeapNode
{
	double m_label;
	int m_index;

public:
	HeapNode() : FibHeapNode() { m_label = 0; };

	virtual void operator =(FibHeapNode& RHS);
	virtual int operator ==(FibHeapNode& RHS);
	virtual int operator <(FibHeapNode& RHS);

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

int HeapNode::operator ==(FibHeapNode& RHS)
{
	if (FHN_compare(RHS)) return 0;
	return m_label == ((HeapNode&) RHS).m_label ? 1 : 0;
}

int HeapNode::operator <(FibHeapNode& RHS)
{
	int X;

	if ((X=FHN_compare(RHS)) != 0) return X < 0 ? 1 : 0;
	return m_label < ((HeapNode&) RHS).m_label ? 1 : 0;
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
	setValues(_D, D, 0);

	int *out_links = (int *) R_alloc(d1*d1, sizeof(int));
	int *out_degree = (int *) R_alloc(d1, sizeof(int));

	for (int i = 0; i < d1; i++)
	{
		out_degree[i] = 0;
		for (int j = 0; j < d1; j++)
		{
			if (W[i + j * d1] != 0)
			{
				out_links[i*d1 + out_degree[i]] = j;
				out_degree[i]++;
			}
		}
	}

	bool *label_free;
	HeapNode *label;

	label_free = (bool *) R_alloc(d1, sizeof(bool));
	label = new HeapNode [d1];

	HeapNode *min_label, temp_label;
	int min_index;
	FibHeap heap;
	double min_value, d;

	for (int s = 0; s < d1; s++)
	{
		// Dijkstra's algorithm

		// initial the nodes
		for (int i = 0; i < d1; i++)
		{
			label[i].setKeyValue(i != s ? HUGE_VAL : 0);
			label[i].setIndexValue(i);
			heap.insert(&label[i]);
			label_free[i] = true;
		}

		// do while
		while (true)
		{	
			min_label = (HeapNode *) heap.extractMin();				
			if (min_label == NULL || min_label->getKeyValue() >= HUGE_VAL)
				break;
			min_index = min_label->getIndexValue();
			D[s + min_index * d1] = min_value = min_label->getKeyValue();
			label_free[min_index] = false;
				
			for (int j = 0; j < out_degree[min_index]; j++)
			{
				int k = out_links[min_index * d1 + j];
				if (label_free[k]) // i, k is neighboring
				{
					d = min_value + W[min_index + k * d1];
					if (d < label[k].getKeyValue())
					{
						temp_label = label[k];
						temp_label.setKeyValue(d);
						heap.decreaseKey(&label[k], temp_label);
					}
				}
			}
		}
	}

	delete[] label;

	UNPROTECT(2);
	return (_D);
}