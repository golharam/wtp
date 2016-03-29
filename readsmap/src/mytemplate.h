#ifndef __mytemplate_h__
#define __mytemplate_h__

//In a list, remove items that precede another in front of it in the list
template <class FwdIt, class pred>
int Inplace_compress(FwdIt first, FwdIt last, pred precede) 
{
    FwdIt a;
    FwdIt c;
    for (c = first, a =c+1; a < last; a++) {
	if (!precede(*c, *a)) continue;
	c++;
	*c = *a;
    }
    return c-first+1;
}

// For two nodes A and B, if A <1 B and B <2 A, then we say A dominate B;
// This algorithm remove all the nodes that are dominated by another.
// It return a list of the nodes left in order of <1.
// By definition the list is also in order of <2.

template <class FwdIt, class pred1, class pred2>
int sort_compress(FwdIt first, FwdIt last, pred1 p1, pred2 p2)
{
    sort(first, last, p1);
    return Inplace_compress(first, last, p2);
}

//
template <class FwdIt1, class FwdIt2, class OutIt, class pred1, class pred2>
int dominate_merge(FwdIt1 first1, FwdIt1 last1, FwdIt2 first2, FwdIt2 last2, OutIt out, pred1 p1, pred2 p2)
{
    OutIt a = merge(first1, last1, first2, last2, out, p1);
    int i = (last1-first1)+(last2-first2);
    return Inplace_compress(out, out+i, p2);
}

//assuming a and b are linked list that already sorted by pred, merge them
// together and return a pointer to the first element
template <class T, class pred>
T linked_merge(T a, T b, pred precede) 
{
    T r, t, w;
    if (a == NULL) return b;
    if (b == NULL) return a;
    r = t = NULL;
    while (a && b) {
	if (precede(a, b)) {
	    if (t) {
		t->set_next(a);
	    } else {
		r = a;
	    }
	    t = a;
	    a = a->get_next();
	} else {
	    if (t) {
		t->set_next(b);
	    } else {
		r = b;
	    }
	    t = b;
	    b = b->get_next();
	}
    }
    if (a) t->set_next(a);
    else t->set_next(b);
    return r;
}

template <class T>
void linkeddisplay(T a)
{
    T b;
    for (b = a; b; b = b->get_next()) {
	b->display();
    }
}

template <class T>
void linkeddisplay2(T first, T last)
{
    T b;
    for (b = first; b; b = b->get_next()) {
	b->display();
	if (b == last) break;
    }
}

template <class T>
int checkforloop(T a)
{
    T b, c;
    for (b = c = a; b; ) {
	b = b->get_next();
	if (b) b = b->get_next();
	c = c->get_next();
	if (c == b) {
	    printf("loop detected\n");
	    linkeddisplay2(c->get_next(), b);
	    return 1;
	}
    }
    return 0;
}

template <class T>
int preddec(const T  a, const T b)
{
    return (a->index() > b->index());
}

template <class T>
int predinc(const T  a, const T b)
{
    return (a->index() < b->index());
}

// This hash_table need same(T *a), and hashvalue() plus a public .next field
// for the class T.
template <class T>
class _hash_table {
 public:
    _hash_table(int t) {
	tablesize = t;
	table = new T*[tablesize];
	CHECK_MEM(table);
        memset(table, 0, sizeof(T *)*tablesize);
    }
    ~_hash_table() {
	CHECK_DEL_ARRAY(table);
    }
    void refresh() {
	memset(table, 0, sizeof(T *)*tablesize);
    }
    T *find(int hashvalue) {
	return table[hashvalue];
    }
    T *find(T *a) {
	int i = a->hashvalue(tablesize);
	T *b;
	for (b = table[i]; b; b= b->next) {
	    if (b->same(a)) return b;
	}
	return NULL;
    }
    int is_in(T *a) {
	T *b = find(a);
	if (b) return 1;
	return 0;
    }
    void insert(T *a) {
	int i = a->hashvalue(tablesize);
	a->next = table[i];
	table[i] = a;
    }
    void insert(T *a, int hash) {
	a->next = table[hash];
	table[hash] = a;
    }
    int delete_one(T *a) { // a is supposed to be in the table
	int i = a->hashvalue(tablesize);
	if (table[i] == a) {
	    table[i] = a->next;
	    return 1;
	}
	T *b = table[i];
	while (b->next && b->next != a) {
	    b = b->next;
	}
	if (b->next) {
	    b->next = a->next;
	    return 1;
	}
	return 0;
    }
    T *return_list() {
	T *first = NULL;
	T *a;
	int i;
	for (i = 0; i < tablesize; i++) {
	    if (table[i]) {
		first = table[i];
		break;
	    }
	}
	for (i++; i < tablesize; i++) {
	    if (table[i] == NULL) continue;
	    for (a = table[i]; a->next; a = a->next) 
		;
	    a->next = first;
	    first = table[i]; 
	}
	return first;
    }
    void delete_itemlist(T *p) {
	while (p) {
	    T *q = p->next;
	    delete p;
	    p = q;
	}
    }
    void delete_items() {
	delete_itemlist(return_list());
    }
 protected:
    T **table;
    int tablesize;
};

template <class T>
T *lastoflist(T *p)
{
    if (!p) return NULL;
    while (1) {
        T *q = p->next;
        if (q) p = q;
        else return p;
    }
}

template <class  T>
void freelist(T *p)
{
    while (p) {
        T *q = p->next;
        delete p;
        p = q;
    }
}

template <class T>
int num_of(T *p)
{
    int i = 0;
    while (p) {
	i++;
	p = p->next;
    }
    return i;
}
#endif
