#ifndef _QUBIC_H
#define _QUBIC_H

#include "struct.h"

/* must be defined */
bits16 **profile;
int col_width;
Edge **edge_list;
Edge *edge_ptr;
char *SY_GETLINE = NULL;
double VER = 1.0;
static int bb[USHRT_MAX];
static char *atom = NULL;
static char delims[] = "\t\r\n";
#define MAXC 100000
static const char USAGE[] = 
"\n===================================================================\n\
[Usage]\n\
qubic(data, [argument list]);\n\
like :\n\
qubic(data, file = 'rQUBIC', q = 0.06, c = 0.95, f = 1, k = 2, r = 1, o = 100, d = 'F')\n\
===================================================================\n\
[Input]\n\
-file : input file must be one of two tab-delimited formats\n\
  A) continuous data (default, use pre-set discretization (see -q and -r))\n\
     -------------------------------------\n\
     o        cond1    cond2    cond3\n\
     gene1      2.4      3.5     -2.4\n\
     gene2     -2.1      0.0      1.2\n\
     -------------------------------------\n\
  B) discrete data with arbitray classes (turn on -d)\n\
     use '0' for missing or insignificant data\n\
     -------------------------------------\n\
     o        cond1    cond2    cond3\n\
     gene1        1        2        2\n\
     gene2       -1        2        0\n\
     -------------------------------------\n\
-q : use quantile discretization for continuous data\n\
     default: 0.06 (see details in Method section in paper)\n\
-r : the number of ranks as which we treat the up(down)-regulated value\n\
     when discretization\n\
     default: 1\n\
-d : discrete data, where user should send their processed data\n\
     to different value classes, see above\n\
-C : the flag using the lower bound of condition number (5 persents of the gene number)\n\
===================================================================\n\
[Output]\n\
-o : number of blocks to report, default: 100\n\
-f : filtering overlapping blocks,\n\
     default: 1 (do not remove any blocks)\n\
-k : minimum column width of the block,\n\
     default: 5% of columns, minimum 2 columns\n\
-c : consistency level of the block (0.5-1.0], the minimum ratio between the\n\
     number of identical valid symbols in a column and the total number \n\
     of rows in the output\n\
     default: 0.95\n\
===================================================================\n";
int r_puts()
{
  puts(USAGE);
  return 1;
}
continuous** alloc2d(int rr, int cc)
{
  continuous** result;
	int i;
	AllocArray(result, rr);
	for (i = 0; i < rr; i++)
		AllocArray(result[i], cc);
	return result;
}
discrete** alloc2c(int rr, int cc)
{
	discrete** result;
	int i;
	AllocArray(result, rr);
	for (i = 0; i < rr; i++)
		AllocArray(result[i], cc);
	return result;
}
char** alloc2c_c(int rr, int cc)
{
  char** result;
	int i;
	AllocArray(result, rr);
	for (i = 0; i < rr; i++)
		AllocArray(result[i], cc);
	return result;
}
static int charset_add(discrete *ar, discrete s)
{
	/*A signed short can hold all the values between SHRT_MIN  and SHRT_MAX inclusive.SHRT_MIN is required to be -32767 or less,SHRT_MAX must be at least 32767*/
	int ps = s + SHRT_MAX;
	if (bb[ps]<0)
	{
		bb[ps] = sigma; 
		ar[sigma++] = s;
	}
	return bb[ps];
}
static continuous quantile_from_sorted_data(const continuous sorted_data[], size_t n, double f)
{
	/*floor function returns the largest integral value less than or equal to x*/
	int i = floor((n-1)*f);
	continuous delta = (n-1)*f-i;
	return (1-delta)*sorted_data[i]+delta*sorted_data[i+1];
}
discrete dis_value(float current, int divided, float *small, int cntl, float *big, int cntu)
{
	int i;
	float d_space = 1.0 / divided;	
	for(i=0; i < divided; i++)
	{		
            if ((cntl > 0) && (current <= quantile_from_sorted_data(small, cntl, d_space * (i+1)))) 
		    return -i-1;
            if ((cntu > 0) && (current >= quantile_from_sorted_data(big, cntu, 1.0 - d_space * (i+1)))) 
		    return i+1;
	}
	return 0;
}
static int compare_continuous (const void *a, const void *b)
{
    const continuous *da = a;
    const continuous *db = b;
    /*make qsort in the increasing order*/
    return (*da < *db)?-1:(*da != *db);
}
void discretize (const char* stream_nm)
{
  int row, col;
	continuous rowdata[cols];
	float big[cols], small[cols];
	int i,cntu,cntl;
	float f1,f2,f3, upper, lower;
	FILE *fw;
	fw = mustOpen(stream_nm, "w");
	for (row = 0; row < rows; row++)
	{
		for (col = 0; col < cols; col++) 
			rowdata[col] = arr[row][col];
		qsort(rowdata, cols, sizeof *rowdata, compare_continuous);
		f1 = quantile_from_sorted_data(rowdata,cols,1-po->QUANTILE); 
		f2 = quantile_from_sorted_data(rowdata,cols,po->QUANTILE);
		f3 = quantile_from_sorted_data(rowdata,cols,0.5);
		if ((f1-f3)>=(f3-f2))
		{
			upper = 2*f3-f2; 
			lower = f2;
		}
		else
		{
			upper = f1; 
			lower = 2*f3-f1;
		}
		cntu = 0; cntl = 0;
		for (i = 0; i < cols; i++)
		{
      if (rowdata[i] < lower) 
			{ 
				small[cntl] = rowdata[i]; 
				cntl++; 
			}
			if (rowdata[i] > upper) 
			{ 
				big[cntu] = rowdata[i]; 
				cntu++; 
			}
		}
		for (col = 0; col < cols; col++)
			arr_c[row][col] = charset_add(symbols, dis_value(arr[row][col],po->DIVIDED, small, cntl, big, cntu));
		fprintf(fw, "row %s :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", genes[row], lower, upper,cntl, cntu);
    /*
    rrules[row + 0*row] = lower;
    rrules[row + 1*row] = upper;
    rrules[row + 2*row] = cntl;
    rrules[row + 3*row] = cntu;
    */
  }
  progress ("Discretization rules are written to %s", stream_nm);
	fclose(fw);
}

void write_imported (const char* stream_nm)
{
	int row, col;
	FILE *fw;
	fw = mustOpen(stream_nm, "w"); 
	fprintf(fw, "o");
	for (col = 0; col < cols; col++)
		fprintf(fw,"\t%s",conds[col]);
	fputc('\n',fw);
	for (row = 0; row < rows; row++)
	{
		fprintf(fw, "%s", genes[row]);
		for (col = 0; col < cols; col++)
      fprintf(fw, "\t%d", symbols[arr_c[row][col]]);
		fputc('\n', fw);
	}
	progress("Formatted data are written to %s", stream_nm);
	fclose(fw);
}
void read_list (FILE* fp)
{
	int i=0, j=0;
	sub_genes_row = 0;
	char line[MAXC];
	while(fgets(line,MAXC,fp)!=NULL)
	{
		atom = strtok(line, delims);
		strcpy(sub_genes[sub_genes_row], atom);	
		sub_genes_row++;
	}

	/*update the sub_list*/
	AllocArray(sublist,rows);
	for (i = 0; i<rows; i++)
		sublist[i] = FALSE;
	for (i=0; i<sub_genes_row; i++)
		for (j=0; j<rows; j++)
			if (strcmp (sub_genes[i], genes[j])==0)
				sublist[j] = TRUE;
}

/*add write_block subroutine prototypes */
/*
char blockfn[LABEL_LEN];
FILE *fw;
*/

void seed_update (const discrete *s)
{
	int i;
	for (i = 0; i < cols; i++) 
		profile[i][s[i]]++;
}

/* scan through all columns and identify the set within threshold,
 * "fuzziness" of the block is controlled by TOLERANCE (-c)
 */
void scan_block (struct dyStack *gene_set, Block *b_ptr)
{
	int i, j;
	int block_rows, cur_rows;
	block_rows = cur_rows = dsSize(gene_set);
	
	int k;
	for (j = 0; j < cols; j++)
		for (k=0; k<sigma; k++) 
			profile[j][k] = 0;
	for (j = 0; j< cur_rows ; j++)
		seed_update(arr_c[dsItem(gene_set,j)]);

	int btolerance = ceil(po->TOLERANCE* block_rows);
	for (j = 0; j < cols; j++)
	{
		/* See if this column satisfies tolerance */
		/* here i start from 1 because symbols[0]=0 */
		for (i = 1; i < sigma; i++)
		{
			if ((profile[j][i] >= btolerance))
			{
				dsPush(b_ptr->conds, j); break;
			}
		}		
	}
	b_ptr->block_cols = dsSize(b_ptr->conds);
}

/*************************************************************************/

/* Identified clusters are backtraced to the original data, by
 * putting the clustered vectors together, identify common column
 */
void print_bc (FILE* fw, Block* b, int num)
{	
	int i, j;
	int block_rows, block_cols;
	int num_1=0,num_2=0;	
	/* block height (genes) */
	block_rows = b->block_rows;
	/* block_width (conditions) */
	block_cols = b->block_cols;
	fprintf(fw, "BC%03d\tS=%d\tPvalue:%LG \n", num, block_rows * block_cols, b->pvalue);
	/* fprintf(fw, "BC%03d\tS=%d\tPvalue:%lf \n", num, block_rows * block_cols, (double)b->pvalue); */
  fprintf(fw, " Genes [%d]: ", block_rows);
	for (i=0; i<dsSize(b->genes); i++)
		fprintf(fw, "%s ", genes[dsItem(b->genes, i)]);
	fprintf(fw, "\n");

	fprintf(fw, " Conds [%d]: ", block_cols);
	for (i=0; i<dsSize(b->conds); i++)
		fprintf(fw, "%s ", conds[dsItem(b->conds, i)]);
	fprintf(fw, "\n");	
	/* the complete block data output */
	for (i=0; i<dsSize(b->genes); i++)
	{
		fprintf(fw,"%10s:",genes[dsItem(b->genes, i)]);
		for (j=0; j<dsSize(b->conds); j++)
		{
			fprintf(fw, "\t%d", symbols[arr_c[dsItem(b->genes, i)][dsItem(b->conds, j)]]);
			if (i==0)
			{
				if (symbols[arr_c[dsItem(b->genes, i)][dsItem(b->conds, j)]] == 1) num_1++;
				if (symbols[arr_c[dsItem(b->genes, i)][dsItem(b->conds, j)]] == -1) num_2++;
			}
		}
		fputc('\n', fw);
		if (i == b->block_rows_pre -1) 
			fputc('\n',fw);
	}
	/*printf ("BC%03d: #of 1 and -1 are:\t%d\t%d\n",num,num_1,num_2);
	fputc('\n', fw);*/
}


/*add cluster subroutine prototypes */
/*
struct rowMatch
{
    int row_id;
    int matches;
};*/

/* Initialize seed */
static int compare_int (const void *a, const void *b)
{
	const int *da = a;
	const int *db = b;    
	return (*da < *db)?-1:(*da != *db);
}
static void update_colcand(bool *colcand, discrete *g1, discrete *g2)
{
	int i;
	for (i=0; i< cols; i++)
		if (colcand[i] && (g1[i] != g2[i])) 
			colcand[i] = FALSE;
}
static int intersect_row(const bool *colcand, discrete *g1, discrete *g2)
/*caculate the weight of the edge with two vertices g1 and g2*/
{
	int i;
	int cnt = 0;
	for (i=0; i< cols; i++)
		if (colcand[i] && (g1[i] == g2[i]) && (g1[i]!=0)) 
			cnt++;
	return cnt;
}
static int reverse_row(const bool *colcand, discrete *g1, discrete *g2)
{
/*caculate the negative correlation between g1 and g2*/	
	int i;
	int cnt = 0;
	for (i = 0; i < cols; i++)
	{
		if (colcand[i] && (symbols[g1[i]] == -symbols[g2[i]])) 
			cnt++;
	}
	return cnt;
} 
static void seed_current_modify (const discrete *s, bool *colcand, int* cnt, int components)
/* calculate the coverage of any row to the current consensus
 * cnt = # of valid consensus columns
 */
{
	int i, k, flag, n;
	int threshold = ceil(components * po->TOLERANCE);
	discrete ss;
	*cnt = 0;
	for (i=0; i<cols; i++)
	{
		flag = 0; ss = s[i];
		for (k=1; k<sigma; k++)
		{
			n = profile[i][k];
			if (k == ss) 
				n++;
			if (n >= threshold)
			{
				flag = k; 
				break;
			}
		}
		if (flag)	
		{
			(*cnt)++;
			colcand[i] = TRUE;
		}
	}
}

static bool check_seed(Edge *e, Block **bb, const int block_id)
/*check whether current edge can be treat as a seed*/
{
	int profiles[rows];
	int i,b1,b2,b3;
	bool fg = FALSE;
	b1 = b2 = -1;
	for (i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes,e->gene_one) && isInStack(bb[i]->genes, e->gene_two) ) 
			return FALSE; 

	for ( i = 0; i < rows; i++) profiles[i] = 0;
	fg = FALSE;	
	for ( i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes, e->gene_one) ) 
		{ 
			fg = TRUE;
		       	break; 
		}
	if (fg) 
		b1 = i;
	fg = FALSE;	
	for ( i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes, e->gene_two) ) 
		{ 
			fg = TRUE; 
			break; 
		}
	if (fg) 
		b2 = i;
	if ( (b1 == -1)||(b2 == -1) ) 
		return TRUE;
	else
	{
		for ( i = 0; i < bb[b1]->block_rows; i++)
			profiles[dsItem(bb[b1]->genes,i)]++;
		for ( i = 0; i < bb[b2]->block_rows; i++)
			profiles[dsItem(bb[b2]->genes,i)]++;
		for ( i = 0; i < rows; i++)
 			if (profiles[i] > 1) 
				return FALSE;
		b3 = MAX(bb[b1]->block_cols, bb[b2]->block_cols);
		if ( e->score <b3/* (bb[b1]->block_cols + bb[b2]->block_cols) / 2*/ ) 
			return FALSE;
		else 
			return TRUE;
	}
	err("never see this message\n");
	return FALSE;
}

long double get_pvalue (continuous a, int b)
{
	int i =0;
	long double one = 1, pvalue=0;
	long double poisson=one/exp(a);
	for (i=0;i<b+300;i++)
	{
		if (i>(b-1)) 
			pvalue=pvalue+poisson;
		else 
			poisson=poisson*a/(i+1);
	}
	return pvalue;
}

static void block_init(Edge *e, Block *b, 
                     struct dyStack *genes, struct dyStack *scores,
                     bool *candidates, const int cand_threshold,
                     int *components, struct dyStack *allincluster, long double *pvalues)
{
	int i,score,top;
	int cnt = 0, cnt_all=0, pid=0;
	continuous cnt_ave=0, row_all = rows;
	long double pvalue;
	int max_cnt, max_i;
	int *arr_rows, *arr_rows_b;
	AllocArray(arr_rows, rows);
	AllocArray(arr_rows_b, rows);	
	bool *colcand;
	AllocArray(colcand, cols);
	for (i=0; i< cols; i++) 
		colcand[i] = FALSE;
	discrete *g1, *g2;
	g1 = arr_c[dsItem(genes,0)];
	g2 = arr_c[dsItem(genes,1)];
	for (i=0; i< cols; i++)
		if ((g1[i] == g2[i])&&(g1[i]!=0)) 
			colcand[i] = TRUE;

	for (i = 0; i < rows; i++)
	{
		arr_rows[i] = intersect_row(colcand, arr_c[dsItem(genes,0)], arr_c[i]);
		arr_rows_b[i] = arr_rows[i];
	}
	/*we just get the largest 100 rows when we initial a bicluster because we believe that 
	 * the 100 rows can characterize the structure of the bicluster 
	 * btw, it can reduce the time complexity*/
	if (rows > 100)
	{		
		qsort(arr_rows_b, rows, sizeof *arr_rows, compare_int);
		top = arr_rows_b[rows -100];
		for (i = 0; i < rows; i++)
			if (arr_rows[i] < top) 
				candidates[i] = FALSE;
	}
	/*calculate the condition low bound for current seed*/
	int cutoff = floor (0.05*rows);
	b->cond_low_bound = arr_rows_b[rows-cutoff];


	while (*components < rows)
	{
		max_cnt = -1;
		max_i = -1;
		(*components)++;
		cnt_all =0;
		cnt_ave = 0;
		/******************************************************/
		/*add a function of controling the bicluster by pvalue*/
		/******************************************************/
		for (i=0; i< rows; i++)
		{
			if (!candidates[i]) continue;
			if (po->IS_list && !sublist[i]) continue;
			cnt = intersect_row(colcand,arr_c[dsItem(genes,0)],arr_c[i]);
			cnt_all += cnt;
			if (cnt < cand_threshold) 
				candidates[i] = FALSE;
			if (cnt > max_cnt)
			{
				max_cnt = cnt;
				max_i = i;
			}
		}
		cnt_ave = cnt_all/row_all;
		pvalue = get_pvalue (cnt_ave, max_cnt);
		if (po->IS_cond)
		{
			if (max_cnt < po->COL_WIDTH || max_i < 0|| max_cnt < b->cond_low_bound) break;
		}
		else
		{
			if (max_cnt < po->COL_WIDTH || max_i < 0) break;
		}
		if (po->IS_area)
			score = *components*max_cnt;
		else
			score = MIN(*components, max_cnt);
		if (score > b->score) 
			b->score = score;
		if (pvalue < b->pvalue)
			b->pvalue = pvalue;
		dsPush(genes, max_i);
		dsPush(scores,score);
		pvalues[pid++] = pvalue;
		update_colcand(colcand,arr_c[dsItem(genes,0)], arr_c[max_i]);
		candidates[max_i] = FALSE;
	}
	/*be sure to free a pointer when you finish using it*/
	free(colcand);
}

static void print_params(FILE *fw)
{
	char filedesc[LABEL_LEN];
	strcpy(filedesc, "continuous");
	if (po->IS_DISCRETE) 
		strcpy(filedesc, "discrete");
	fprintf(fw, "# QUBIC version %.1f output\n", VER);
	fprintf(fw, "# Datafile %s: %s type\n", po->FN, filedesc);
	fprintf(fw, "# Parameters: -k %d -f %.2f -c %.2f -o %d",
			po->COL_WIDTH, po->FILTER, po->TOLERANCE, po->RPT_BLOCK);
	if (!po->IS_DISCRETE) 
		fprintf(fw, " -q %.2f -r %d", po->QUANTILE, po->DIVIDED);
	fprintf(fw, "\n\n");
}

static int block_cmpr(const void *a, const void *b)
/* compare function for qsort, descending by score */
{
	return ((*(Block **)b)->score - (*(Block **)a)->score);
}

static void sort_block_list(Block **el, int n)
{
	qsort(el, n, sizeof *el, block_cmpr);
}

/************************************************************************/
static int report_blocks(FILE* fw, Block** bb, int num)
{
	print_params(fw);
	sort_block_list(bb, num);
	
	int i, j,k;
	/*MIN MAX et al functions can be accessed in struct.h*/
        int n = MIN(num, po->RPT_BLOCK);
	bool flag;

	Block **output;
	AllocArray(output, n);

	Block **bb_ptr = output;
	Block *b_ptr;
	double cur_rows, cur_cols;
	double inter_rows, inter_cols;
        /*double proportion;*/
	
	/* the major post-processing here, filter overlapping blocks*/
	i = 0; j = 0;
	while (i < num && j < n)
	{
		b_ptr = bb[i];
		cur_rows = b_ptr->block_rows;
		cur_cols = b_ptr->block_cols;

		flag = TRUE;
		k = 0;
		while (k < j)
		{
			inter_rows = dsIntersect(output[k]->genes, b_ptr->genes);
			inter_cols = dsIntersect(output[k]->conds, b_ptr->conds);
			
			if (inter_rows*inter_cols > po->FILTER*cur_rows*cur_cols)
			{
				flag = FALSE; 
				break;
			}
                        k++;
		}
	        i++;
		if (flag)
		{
			print_bc(fw, b_ptr, j++);
			*bb_ptr++ = b_ptr;
		}
	}
	return j;
}

/************************************************************************/

int cluster (FILE *fw, Edge **el, int n)
{
	int block_id = 0;
	Block **bb;
	int allocated = po->SCH_BLOCK;
	AllocArray(bb, allocated);

	Edge *e;
	Block *b;
	struct dyStack *genes, *scores, *b_genes, *allincluster;
	
	int i, j, k, components;
	AllocArray(profile, cols);
	for (j = 0; j < cols; j++) 
		AllocArray(profile[j], sigma);

	genes = dsNew(rows);
	scores = dsNew(rows);
	allincluster = dsNew(rows);

    
	long double *pvalues;
	AllocArray(pvalues, rows);

	bool *candidates;
	AllocArray(candidates, rows);

	e = *el; 
	i = 0;
	while (i++ < n)
	{	
		e = *el++;
		/* check if both genes already enumerated in previous blocks */
		bool flag = TRUE;
		/* speed up the program if the rows bigger than 200 */
	        if (rows > 250)
		{ 
			if ( isInStack(allincluster,e->gene_one) && isInStack(allincluster,e->gene_two) )
				flag = FALSE;
			else if ((po->IS_TFname)&&(e->gene_one!= TFindex)&&(e->gene_two!=TFindex))
				flag = FALSE;
			else if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two]))
				flag =FALSE;
		}
		else   
		{
			flag = check_seed(e, bb, block_id);
			if ((po->IS_TFname)&&(e->gene_one!= TFindex)&&(e->gene_two!=TFindex))
				flag = FALSE;
			if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two]))
				flag = FALSE;
		}
		if (!flag) continue;

		for (j = 0; j < cols; j++)
			for (k = 0; k < sigma; k++) 
				profile[j][k] = 0;

		/*you must allocate a struct if you want to use the pointers related to it*/
		AllocVar(b);
		/*initial the b->score*/
                b->score = MIN(2, e->score);
		/*initial the b->pvalue*/
		b->pvalue = 1;
	
		/* initialize the stacks genes and scores */		
		int ii;		
		dsClear(genes);
		dsClear(scores);		
		for(ii = 0; ii < rows; ii ++)
		{
			dsPush(genes,-1);
			dsPush(scores,-1);
		}		
		dsClear(genes);
		dsClear(scores);
		
		dsPush(genes, e->gene_one);
		dsPush(genes, e->gene_two);
		dsPush(scores, 1);
		dsPush(scores, b->score);

		/* branch-and-cut condition for seed expansion */
		int cand_threshold = floor(po->COL_WIDTH * po->TOLERANCE);
                if (cand_threshold < 2) 
			cand_threshold = 2;

		/* maintain a candidate list to avoid looping through all rows */		
		for (j = 0; j < rows; j++) 
			candidates[j] = TRUE;
		candidates[e->gene_one] = candidates[e->gene_two] = FALSE;
		components = 2;

		/* expansion step, generate a bicluster without noise */
		block_init(e, b, genes, scores, candidates, cand_threshold, &components, allincluster, pvalues);

		/* track back to find the genes by which we get the best score*/
		for(k = 0; k < components; k++)
		{
			if (po->IS_pvalue)
				if ((pvalues[k] == b->pvalue) &&(k >= 2) &&(dsItem(scores,k)!=dsItem(scores,k+1))) break;
			if ((dsItem(scores,k) == b->score)&&(dsItem(scores,k+1)!= b->score)) break;
		}
		components = k + 1;
		int ki;
		for (ki=0; ki < rows; ki++)
			candidates[ki] = TRUE;

		for (ki=0; ki < components - 1 ; ki++)
		{
			seed_update(arr_c[dsItem(genes,ki)]);
			candidates[dsItem(genes,ki)] = FALSE;
		}
		candidates[dsItem(genes,k)] = FALSE;
		genes->top = k ;
		int cnt = 0;
		bool *colcand;
		AllocArray(colcand, cols);
		for(ki = 0; ki < cols; ki++) 
			colcand[ki] = FALSE;             
    
		/* add columns satisfy the conservative r */ 
		seed_current_modify(arr_c[dsItem(genes,k)], colcand, &cnt, components);
		
		/* add some new possible genes */
		int m_cnt;
		for ( ki = 0; ki < rows; ki++)
		{
			if (po->IS_list && !sublist[ki]) continue;
			m_cnt = intersect_row(colcand, arr_c[dsItem(genes,0)], arr_c[ki]);
			if ( candidates[ki] && (m_cnt >= floor(cnt* po->TOLERANCE)) )
			{
				dsPush(genes,ki);
				components++;
				candidates[ki] = FALSE;
			}
		}
                b->block_rows_pre = components;
		
		/* add genes that negative regulated to the consensus */
		for ( ki = 0; ki < rows; ki++)
		{
			if (po->IS_list && !sublist[ki]) continue;
			m_cnt = reverse_row(colcand, arr_c[dsItem(genes,0)], arr_c[ki]);
			if ( candidates[ki] && (m_cnt >= floor(cnt * po->TOLERANCE)) )
			{
				dsPush(genes,ki);
				components++;
				candidates[ki] = FALSE;
			}
		}
		free(colcand);

		/* save the current cluster*/
		b_genes = dsNew(b->block_rows_pre);
		for (ki = 0; ki < b->block_rows_pre; ki++)
			dsPush(b_genes, dsItem(genes,ki));

		/* store gene arrays inside block */
		b->genes = dsNew(components);
		b->conds = dsNew(cols);
	
		scan_block(b_genes, b);
		if (b->block_cols == 0) continue;
		b->block_rows = components;
                if (po->IS_pvalue)
			b->score = -(100*log(b->pvalue));
		else
			b->score = b->block_rows * b->block_cols;		

		dsClear(b->genes);
		for ( ki=0; ki < components; ki++)
			dsPush(b->genes,dsItem(genes,ki));
		for(ki = 0; ki < components; ki++)
			if(!isInStack(allincluster, dsItem(genes,ki))) 
				dsPush(allincluster,dsItem(genes,ki));	
		/*save the current block b to the block list bb so that we can sort the blocks by their score*/
		bb[block_id++] = b;

		/* reaching the results number limit */
		if (block_id == po->SCH_BLOCK) break;
		verboseDot();	
	}
	/* writes character to the current position in the standard output (stdout) and advances the internal file position indicator to the next position.
	 * It is equivalent to putc(character,stdout).*/
	putchar('\n');
	/* free-up the candidate list */
	free(candidates);
	free(allincluster);
	free (pvalues);
	return report_blocks(fw, bb, block_id);
}

/*make_graph subroutine prototypes */
// static const int HEAP_SIZE = 20000000;
static const int HEAP_SIZE = 20000000;

void seed_deduct (const discrete *s)
/* remove a row from the profile */
{
	int i;
	discrete ss;
	for (i=0; i<cols; i++)
	{
		ss = s[i];
		profile[i][ss]--;
	}
}

static int str_intersect_r (const discrete *s1, const discrete *s2)
{
  int common_cnt = 0;
	/* s1 and s2 of equal length, so we check s1 only */
	int i;
	for (i=0;i<cols;i++)
		if (s1[i]==s2[i] && (s1[i]!=0)) 
			common_cnt++;
	return common_cnt;
}

static int edge_cmpr(void *a, void *b)
{
	int score_a, score_b;
	score_a = ((Edge *)a)->score;
	score_b = ((Edge *)b)->score;

	if (score_a < score_b) return -1;
	if (score_a == score_b) return 0;
	return 1;
}

static void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min)
{
	if (a->fh_n < HEAP_SIZE) 
	{
		fh_insert(a, (void *)i);
	}
	else
	{
		if (edge_cmpr(cur_min, i) < 0)
		{
			/* Remove least value and renew */
			fh_extractmin(a);
			fh_insert(a, (void *)i);
			/* Keep a memory of the current min */
			*cur_min = (Edge *)fh_min(a);
		}
	}
}

static void fh_dump(struct fibheap *a, Edge **res)
{
	int i;
	int n = a->fh_n;
	for (i=n-1; i>=0; i--)
		res[i] = (Edge *) fh_extractmin(a);
}

void make_graph (const char* fn)
{
	FILE *fw = mustOpen(fn, "w");
	int i, j, cnt;
	int rec_num = 0;
	if (po->COL_WIDTH == 2) 
		po->COL_WIDTH = MAX(cols/20, 2);
	
	/* edge_ptr describe edges */
	AllocArray(edge_list, HEAP_SIZE);

	/* Allocating heap structure */
	struct fibheap *heap;
	heap = fh_makeheap();
	fh_setcmp(heap, edge_cmpr);

	/* Generating seed list and push into heap */
	progress("Generating seed list (minimum weight %d)", po->COL_WIDTH);
	Edge __cur_min = {0, 0, po->COL_WIDTH};
	Edge *_cur_min = &__cur_min;
	Edge **cur_min = & _cur_min;
	/* iterate over all genes to retrieve all edges */
	for (i = 0; i < rows; i++)
    for (j = i+1; j < rows; j++)
		{
      cnt = str_intersect_r(arr_c[i], arr_c[j]);
      if (cnt < (*cur_min)->score) continue;
			AllocVar(edge_ptr);
			edge_ptr -> gene_one = i;
			edge_ptr -> gene_two = j;
			edge_ptr -> score = cnt;	
			fh_insert_fixed(heap, edge_ptr, cur_min);
		}
	rec_num = heap->fh_n;
	if (rec_num == 0)
		errAbort("Not enough overlap between genes");
	/* sort the seeds */
	uglyTime("%d seeds generated", rec_num);
	ReAllocArray(edge_list, rec_num);
	fh_dump(heap, edge_list);
	/* bi-clustering */
  int n_blocks = 0;
	progress("Clustering started");
	n_blocks = cluster(fw, edge_list, rec_num);
	uglyTime("%d clusters are written to %s", n_blocks, fn);
	/* clean up */
	for (i=0; i<rec_num; i++)
		free(edge_list[i]);
	free(edge_list);
  fclose(fw);
}

/* expand subroutine prototypes */
discrete** another_arr_c;
char** another_genes;
char** another_conds;
int another_rows;
int another_cols;

static int intersect_rowE(const bool *colcand, discrete *g1, discrete *g2, const int cols)
{
	int i, cnt = 0;
	for(i=0; i< cols; i++)
		if( colcand[i] && (g1[i] == g2[i]) && g1[i] != 0 ) 
			cnt++;
	return cnt;
}

static int reverse_rowE(const bool *colcand, discrete *g1, discrete *g2, const int cols)
{
	int i, cnt = 0;
	for(i=0; i< cols; i++)
		if( colcand[i] && (symbols[g1[i]] == -symbols[g2[i]])) 
			cnt++;
	return cnt;
}

void store_block(Block *b_ptr, struct dyStack *ge, struct dyStack *co)
{
	int row, col;
	row = dsSize(ge);
	col = dsSize(co);
	b_ptr->genes = ge;
	b_ptr->conds = co;
	b_ptr->block_rows = row;
	b_ptr->block_cols = col;
}
static void init_expand()
{
	another_genes = genes;
	another_conds = conds;
	another_arr_c = arr_c;
	another_rows = rows;
	another_cols = cols;
}
/* Read the .block file, get components and colcand */
void read_and_solve_blocks(FILE *fb, const char *fn)
{
	init_expand();
	int col;
	char *line = NULL;
	int bnumber = 0;
	struct dyStack *ge, *co;
	int i, components, m_cnt;
	bool *colcand;
	bool *candidates;
	Block *b;
	AllocVar(b);
	AllocArray(colcand, another_cols);
	AllocArray(candidates, another_rows);
	ge = dsNew(another_rows);
	co = dsNew(another_cols);
	FILE *fo = mustOpen(fn, "w");

	/* main course starts here */
	while (fgets(line, MAXC, fb)!=NULL)
	{
	        /* fast forward to a line that contains BC*/
		/* strncmp compares up to num characters of the C string str1 to those of the C string str2
		 * strncmp ( const char * str1, const char * str2, size_t num )*/
		while (strncmp(line, "BC", 2)!=0) 
		{
			if (fgets(line, MAXC, fb)==NULL) 
			{
				uglyTime("expanded biclusters are written to %s", fn);
				exit(0);
			}
		}
		components = 0;
		col = 0;
		dsClear(ge);
		dsClear(co);
		for (i=0; i< another_cols; i++)
			colcand[i] = FALSE;
		for (i=0; i< another_rows; i++)
			candidates[i] = TRUE;
		/* read genes from block */		
		SY_GETLINE = fgets(line, MAXC, fb);
		atom = strtok(line, delims);
		atom = strtok(NULL, delims);
		while((atom = strtok(NULL, delims)) != NULL)
		{
			/* look up for genes number */
			if (strlen(atom) == 0) continue;			
			for(i=0; i<another_rows; i++)
			{
				if (strcmp(atom ,another_genes[i]) == 0) break;
			}
			candidates[i] = FALSE;			
			dsPush(ge, i);
			components++;
		}
		/* read conditions from block */
		SY_GETLINE = fgets(line, MAXC, fb);
		atom = strtok(line, delims);
		atom = strtok(NULL, delims);
		while((atom = strtok(NULL, delims)) != NULL)
		{
			/*if (strlen(atom) < 5) break;*/			
			if (strlen(atom) == 0) continue;			
			for(i=0; i<another_cols; i++)
				if (strcmp(atom, another_conds[i]) == 0) break;
			colcand[i] = TRUE;
			dsPush(co, i);
			col++;
		}
		
		b->block_rows_pre = components;
		/* add some possible genes */
		for( i = 0; i < another_rows; i++)
		{
			m_cnt = intersect_rowE(colcand, another_arr_c[dsItem(ge,0)], another_arr_c[i], another_cols);
			/*printf ("%d\n",m_cnt);*/
			if( candidates[i] && (m_cnt >= (int)floor( (double)col * po->TOLERANCE)) )
			{
				dsPush(ge,i);
				components++;
				candidates[i] = FALSE;
			}
		}
		/* add genes that negative regulated to the consensus */
		for( i = 0; i < another_rows; i++)
		{
			m_cnt = reverse_rowE(colcand, another_arr_c[dsItem(ge,0)], another_arr_c[i], another_cols);
			if( candidates[i] && (m_cnt >= (int)floor( (double)col * po->TOLERANCE)) )
			{
				dsPush(ge,i);
				components++;
				candidates[i] = FALSE;
			}
		}
		if(dsSize(ge) > 1)
		{		
			store_block(b, ge, co);
			/*another_print_bc(fo, b, bnumber);*/
			print_bc(fo, b, bnumber++);
		}
	}
			printf ("1111\n");
}

#if 0
void another_print_bc(FILE* fw, Block* b, int num)
{	
	int i, j;
	int block_rows, block_cols;
	
	/* block height (genes) */
	block_rows = b->block_rows;
	/* block_width (conditions) */
	block_cols = b->block_cols;

    /*if (block_rows==0 || block_cols==0) return;*/
	fprintf(fw, "BC%03d\tS=%d\n", num, block_rows * block_cols);

	fprintf(fw, "Genes [%d]:", block_rows);
	for (i=0; i<dsSize(b->genes); i++)
		fprintf(fw, "%s ", another_genes[dsItem(b->genes, i)]);
	fprintf(fw, "\n");

	fprintf(fw, "Conds [%d]:", block_cols);
	for (i=0; i<dsSize(b->conds); i++)
		fprintf(fw, "%s ", another_conds[dsItem(b->conds, i)]);
	fprintf(fw, "\n");	
	/* the complete block data output */
	for (i=0; i<dsSize(b->genes); i++)
	{
		fprintf(fw,"%10s:",another_genes[dsItem(b->genes, i)]);
		for (j=0; j<dsSize(b->conds); j++)
		{
            fprintf(fw, "\t%d", symbols[arr_c[dsItem(b->genes,i)]][dsItem(b->cons, j)]);
		}
		fputc('\n', fw);
		
		if (i == b->block_rows_pre -1) fputc('\n',fw); /* add a blank line for negative */
	}
	fputc('\n', fw);
}
void get_chars_size(FILE *fp)
{
	size_t n = 0;
	char *line;
	
	getline(&line, &n ,fp);
	
	another_cols = 0;
	atom = strtok(line, delims);
	while (strtok(NULL, delims) != NULL)/* not to the end of line */
	{
		++another_cols;
	}
	another_rows = 0;
	while (getline(&line, &n, fp) >= 0)
	{
		++another_rows;
	}
	fseek(fp, 0, 0);
    progress("File contains %d rows and %d cols", another_rows, another_cols);
}

void read_array(FILE *fp)
{
	size_t n = 0;
	char *line;
	int row = 0;
	int col = 0;

	/* initialization */	
	another_genes = alloc2c(another_rows, LABEL_LEN);
	another_conds = alloc2c(another_cols, LABEL_LEN);
	
    another_arr_c = alloc2c(another_rows, another_cols+1);
	for ( row = 0; row < another_rows; row++)
		another_arr_c[row][another_cols] = '\0';
	
	/* read condition names */
	getline(&line, &n, fp);
	atom = strtok(line, delims);
	while (atom != NULL)
	{
		strcpy(another_conds[col], atom);
		atom = strtok(NULL, delims);
		if (++col == another_cols) break;
	}
	/* read array .chars */
	row = 0;	
	while (getline(&line, &n, fp) >=0)
	{			
		atom = strtok(line, delims);
		strcpy(another_genes[row], atom);
		atom = strtok(NULL, delims);
		for ( col=0; col < another_cols; col++)
			another_arr_c[row][col] = atom[col];
		if (++row == another_rows) break;
	}
	fclose(fp);
}	
#endif

int init_qubic( double *r_data, char **r_rowsnames, char **r_colsnames, int *r_rows, int *r_cols, char *tfile, double *rq, double *rc, double *rf, int *rk, int *rr, int *ro, int *rd )
{
  int i = 0, j = 0;
	rows = *r_rows;
	cols = *r_cols;
	genes = r_rowsnames;
	conds = r_colsnames;
	arr = alloc2d(rows, cols);
	arr_c = alloc2c(rows, cols);
	discrete temp_dtoi = 0;
	for(i = 0; i < rows; i++)
    	for(j = 0; j < cols; j++)
    		arr[i][j] = r_data[i+j*rows];
	uglyTime(NULL);
	printf("\nQUBIC %.1f: greedy biclustering (compiled "__DATE__" "__TIME__")\n\n", VER);
	/* get the program options defined in get_options.c */
	/*set memory for the point which is decleared in struct.h*/
	AllocVar(po);
	/*Initialize the point*/
	strcpy(po->FN, tfile);
	strcpy(po->BN, " ");
	strcpy(po->LN, " ");
	/* case 'l': strcpy(po->LN, optarg); po->IS_list =TRUE; */
	strcpy(po->TFname, " ");
	/* case 'T': strcpy(po->TFname, optarg); po->IS_TFname = TRUE; */
  if(*rd == 1)  po->IS_DISCRETE = TRUE;
  else  po->IS_DISCRETE = FALSE;
  // po->IS_DISCRETE = TRUE;
	po->IS_TFname = FALSE;
	po->IS_pvalue = FALSE;
	/* case 'P': po->IS_pvalue = TRUE; */
	po->COL_WIDTH = *rk;
	po->DIVIDED = *rr;
	po->QUANTILE = *rq;
	po->TOLERANCE = *rc;
	po->FP = NULL;
	po->FB = NULL;
	po->RPT_BLOCK = *ro;
	po->SCH_BLOCK = 2*po->RPT_BLOCK;
	/* ensure enough searching space */
	/*if (po->SCH_BLOCK < 1000) po->SCH_BLOCK = 1000;*/ 
	po->FILTER = *rf;
	po->IS_SWITCH = FALSE;
	/* case 's': po->IS_SWITCH = TRUE; */
	po->IS_area = FALSE;
	/* case 'S': po->IS_area = TRUE; */
	po->IS_cond = FALSE;
	/* case 'C': po->IS_cond = TRUE; */
	po->IS_list = FALSE;
	if(po->IS_SWITCH) 
	{	
		po->IS_DISCRETE = TRUE; 
		po->FB = mustOpen(po->BN, "r");
	}
	if(po->IS_list)
	{
		po->FL = mustOpen(po->LN, "r");
	}
	/* check if there exist a gene name equals to TFname by -T */
	for(int i = 0; i < rows; i++)
		if(strcmp (genes[i], po->TFname) == 0)
			printf ("%d\n", i);
  AllocArray(symbols, USHRT_MAX);
  memset(bb, -1, USHRT_MAX*sizeof(*bb));
	charset_add(symbols, 0); 
	if(po->IS_DISCRETE)
	{
		for(i = 0; i < rows; i++)
			for(j = 0; j < cols; j++)
			{
				temp_dtoi = (discrete)arr[i][j];
				arr_c[i][j] = charset_add(symbols, temp_dtoi);
			}
		printf("Discretized data contains %d classes with charset [ ", sigma);
		for(i = 0 ; i < sigma; i++)
			printf("%d ", symbols[i]);  printf("]\n");
	}
	else
  {
    for(i = 0; i < rows; i++)
      for(j = 0; j < cols; j++)
        arr_c[i][j] = 0;
    discretize(addSuffix(tfile, ".rules"));
  }
	/*read in the sub-gene list*/
	if(po->IS_list)
	{
		sub_genes = alloc2c_c(rows, LABEL_LEN);
		read_list(po->FL);
	}
	/*we can do expansion by activate po->IS_SWITCH*/
	if(po->IS_SWITCH)  read_and_solve_blocks(po->FB, addSuffix(tfile, ".expansion"));
	else
	{
		/* formatted file */
		write_imported(addSuffix(tfile, ".chars"));
    for (i = 0; i < rows; i++)
		  for (j = 0; j < cols; j++)
        r_data[i+j*rows] = (int)symbols[arr_c[i][j]];
		/* the file that stores all blocks */
		if (po->IS_list)
			make_graph(addSuffix(tfile, ".block"));
		else
			make_graph(addSuffix(tfile, ".blocks"));
 	}
	free(po);
	free(sublist);
	return 1;
}

int cgetbc(double *rbc, int *ro, char **filer, char **filew)
{
  FILE *fpr = fopen(filer[0], "r");
  FILE *fpw = fopen(addSuffix(filew[0], ".bc"), "w");
  int i = 0, n = *ro, num = 0, ntemp = 0, nbc = 0;
  char temp[1000];
  char czero = '0', tempnum[20];
  while(fscanf(fpr, "%s", temp) != EOF)
  {
    if(strcmp(temp, "Genes") == 0)
    {
      fprintf(fpw, "BC%d\n", nbc);
      fscanf(fpr, "%s", tempnum);
      num = 0;
      for(i = 0; i < 20; i++)
      {
        if(tempnum[i] == ']')  break;
        else if(tempnum[i] == ':')  break;
        else if(tempnum[i] == '[')  continue;
        else
        {
          num *= 10;
          ntemp = tempnum[i] - czero;
          num += ntemp;
        }
      }
      /* fprintf(fpw, "%d\n", num); */
      for(i = 0; i < num; i++)
      {
        fscanf(fpr, "%s", temp);
        fprintf(fpw, "%s\t", temp);
      }
      fprintf(fpw, "\n");
      rbc[nbc+n] = num;
      fscanf(fpr, "%s", temp);
      if(strcmp(temp, "Conds") == 0)
      {
        /* fprintf(fpw, "%s\t", temp); */
        fscanf(fpr, "%s", tempnum);
        num = 0;
        for(i = 0; i < 20; i++)
        {
          if(tempnum[i] == ']')  break;
          else if(tempnum[i] == ':')  break;
          else if(tempnum[i] == '[')  continue;
          else
          {
            num *= 10;
            ntemp = tempnum[i] - czero;
            num += ntemp;
          }
        }
        /* fprintf(fpw, "%d\n", num); */
        for(i = 0; i < num; i++)
        {
          fscanf(fpr, "%s", temp);
          fprintf(fpw, "%s\t", temp);
        }
        fprintf(fpw, "\n");
        rbc[nbc+2*n] = num;
        rbc[nbc] = nbc;
        nbc++;
      }
    }
  }
  if(nbc < n)
  {
    nbc = n - nbc;
    for(i = 0; i < nbc; i++)
      rbc[nbc] = -1;
  }
  fclose(fpw);
  fclose(fpr);
  return 1;
}
#endif
