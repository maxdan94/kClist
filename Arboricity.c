/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.
This is an implementation of the 1985 paper of Chiba And Nishizeki "Arboricity and subgraph listing".

To compile:
"gcc Arboricity.c -O9 -o Arboricity".

To execute:
"./Arboricity edgelist.txt k".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
Will print the number of k-cliques.
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>


#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg ;


typedef struct {

	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges

	unsigned *ns;//ns[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *cd;//cumulative degree: (start with 0) length=n+1
	unsigned *adj;//list of neighbors with lower degree

	unsigned char *lab;//lab[i] label of node i
	unsigned **sub;//sub[l]: nodes in G_l

} specialsparse;


void freespecialsparse(specialsparse *g, unsigned char k){
	unsigned char i;
	free(g->ns);
	for (i=2;i<k+1;i++){
		free(g->d[i]);
		free(g->sub[i]);
	}
	free(g->d);
	free(g->sub);
	free(g->cd);
	free(g->adj);
	free(g);
}

//Compute the maximum of three unsigned integers.
inline unsigned int max3(unsigned int a,unsigned int b,unsigned int c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

specialsparse* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	specialsparse *g=malloc(sizeof(specialsparse));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t))==2) {//Add one edge
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		g->e++;
		if (g->e==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;

	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

//Building the special graph structure
void mkspecial(specialsparse *g, unsigned char k){
	unsigned i,ns,max;
	unsigned *d,*sub;
	unsigned char *lab;

	d=calloc(g->n,sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		d[g->edges[i].s]++;
		d[g->edges[i].t]++;
	}

	g->cd=malloc((g->n+1)*sizeof(unsigned));
	ns=0;
	g->cd[0]=0;
	max=0;
	sub=malloc(g->n*sizeof(unsigned));
	lab=malloc(g->n*sizeof(unsigned char));
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		if (d[i-1]>0){
			max=(max>d[i-1])?max:d[i-1];
			sub[ns]=i-1;
			ns++;
			d[i-1]=0;
			lab[i-1]=k;
		}
	}

	g->adj=malloc(2*g->e*sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		g->adj[ g->cd[g->edges[i].s] + d[ g->edges[i].s ]++ ]=g->edges[i].t;
		g->adj[ g->cd[g->edges[i].t] + d[ g->edges[i].t ]++ ]=g->edges[i].s;
	}
	free(g->edges);

	g->ns=malloc((k+1)*sizeof(unsigned));
	g->ns[k]=ns;

	g->d=malloc((k+1)*sizeof(unsigned*));
	g->sub=malloc((k+1)*sizeof(unsigned*));
	for (i=2;i<k;i++){
		g->d[i]=malloc(g->n*sizeof(unsigned));
		g->sub[i]=malloc(max*sizeof(unsigned));
	}
	g->d[k]=d;
	g->sub[k]=sub;

	g->lab=lab;
}

/*For futur use in qsort_r.
int compare_r(void const *p_i, void const *p_j,void *p_deg){
	unsigned i = *((unsigned*)p_i);
	unsigned j = *((unsigned*)p_j);
	unsigned *deg = (unsigned*)p_deg;
	return (deg[i] < deg[j]) ? 1 : -1;
}
*/

void arg_bucket_sort(unsigned *key, unsigned n, unsigned *val) {
  unsigned i,j;
  static unsigned *c=NULL, *cc=NULL, *key2=NULL;
  if (c==NULL){
    c=malloc(n*sizeof(unsigned));//count
    cc=malloc(n*sizeof(unsigned));//cummulative count
    key2=malloc(n*sizeof(unsigned));//sorted array
  }
  bzero(c,n*sizeof(unsigned));

  for (i=0;i<n;i++){
    (c[val[key[i]]])++;
  }
  cc[0]=0;
  for (i=1;i<n;i++){
    cc[i]=cc[i-1]+c[i-1];
    c[i-1]=0;
  }
	c[i-1]=0;

  for (i=0;i<n;i++){
		j=val[key[i]];
    key2[cc[j] + c[j]++]=key[i];
  }

}

void kclique(unsigned l, specialsparse *g, unsigned long long *n) {
	unsigned i,j,k,end,u,v,w;

	if(l==2){
		for(i=0; i<g->ns[2]; i++){//list all edges
			u=g->sub[2][i];
			end=g->cd[u]+g->d[2][u];
			for (j=g->cd[u];j<end;j++) {
				v=g->adj[j];
				if(v<u){
					(*n)++;//listing here!!!
				}
			}
		}
		return;
	}

	arg_bucket_sort(g->sub[l], g->ns[l], g->d[l]);
	//qsort_r(g->sub[l],g->ns[l],sizeof(unsigned),compare_r,g->d[l]);//qsort and bucket sort leads to similar running time in practice

	for(i=0; i<g->ns[l]; i++){
		u=g->sub[l][i];
		g->ns[l-1]=0;
		end=g->cd[u]+g->d[l][u];
		for (j=g->cd[u];j<end;j++){//relabeling nodes and forming U'.
			v=g->adj[j];
			if (g->lab[v]==l){
				g->lab[v]=l-1;
				g->sub[l-1][g->ns[l-1]++]=v;
				g->d[l-1][v]=0;//new degrees
			}
		}
		for (j=0;j<g->ns[l-1];j++){//reodering adjacency list and computing new degrees
			v=g->sub[l-1][j];
			end=g->cd[v]+g->d[l][v];
			for (k=g->cd[v];k<end;k++){
				w=g->adj[k];
				if (g->lab[w]==l-1){
					g->d[l-1][w]++;
				}
				else{
					g->adj[k--]=g->adj[--end];
					g->adj[end]=w;
				}
			}
		}

		kclique(l-1, g, n);

		for (j=0;j<g->ns[l-1];j++){//moving u to last position in each entry of the adjacency list
			v=g->sub[l-1][j];
			g->lab[v]=l;
			end=g->cd[v]+g->d[l-1][v];
			for (k=g->cd[v];k<end;k++){
				w=g->adj[k];
				if(w==u){
					g->adj[k]=g->adj[--end];
					g->adj[end]=w;
					g->d[l-1][v]--;
					break;
				}
			}
		}
		g->lab[u]=l+1;

	}

	return;
}


int main(int argc,char** argv){
	specialsparse* g;
	unsigned char k=atoi(argv[1]);
	unsigned long long n;

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[2]);

	g=readedgelist(argv[2]);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Number of nodes = %u\n",g->n);
	printf("Number of edges = %u\n",g->e);
	printf("Building the graph structure\n");

	mkspecial(g,k);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Iterate over all cliques\n");

	//qsort_r(g->sub[k],g->ns[k],sizeof(unsigned),compare_r,g->d[k]);
	//arg_bucket_sort(g->sub[k], g->ns[k], g->d[k]);
	n=0;
	kclique(k, g, &n);

	printf("Number of %u-cliques: %llu\n",k,n);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	freespecialsparse(g,k);

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}
