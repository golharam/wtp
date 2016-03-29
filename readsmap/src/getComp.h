
class getComp 
{
    public:
	getComp();
	int get_sum(char *seq, int len);
	int is_compatible(char *seq1, int len1, char *seq2, int len2);
    protected:
	char addColors[5][128];
}; 
