/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.
This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing".

To compile:
"gcc kClistEdgeParallel.c -O9 -o kClistEdgeParallel -fopenmp".

To execute:
"./kClistEdgeParallel p k edgelist.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
k is the size of the k-cliques
p is the number of threads
Will print the number of k-cliques.
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <omp.h>


#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed


typedef unsigned long int Node;
typedef unsigned long long int Edge;
typedef unsigned long long int Clique;
typedef unsigned char Kvalue;

typedef struct {
	Node s;
	Node t;
} edge;

typedef struct {
	Node node;
	Node deg;
} nodedeg ;

typedef struct {
	Node n;//number of nodes
	Edge e;//number of edges
	edge *edges;//list of edges
	Node *rank;//ranking of the nodes according to degeneracy ordering
	//unsigned *map;//oldID newID correspondance NOT USED IN THIS VERSION
} edgelist;

typedef struct {
	Node n;
	Edge e;
	edge *edges;//ading this again here: TO IMPROVE
	Edge *cd;//cumulative degree: (starts with 0) length=n+1
	Node *adj;//truncated list of neighbors
	Node core;//core value of the graph
} graph;

typedef struct {
	Node *n;//n[l]: number of nodes in G_l
	Node **d;//d[l]: degrees of G_l
	Node *adj;//truncated list of neighbors
	unsigned char *lab;//lab[i] label of node i
	Node **nodes;//sub[l]: nodes in G_l
	Node core;
	
	Node* new;//forgot what that is...
	Node* old;
} subgraph;

void free_graph(graph *g){
	free(g->cd);
	free(g->adj);
	free(g);
}

void free_subgraph(subgraph *sg, Kvalue k){
	Kvalue i;
	free(sg->n);
	for (i=2;i<k;i++){
		free(sg->d[i]);
		free(sg->nodes[i]);
	}
	free(sg->d);
	free(sg->nodes);
	free(sg->lab);
	free(sg->adj);
	free(sg);
}


//Compute the maximum of three unsigned integers.
inline Node max3(Node a,Node b,Node c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

edgelist* readedgelist(char* input){
	Edge e1=NLINKS;
	edgelist *el=malloc(sizeof(edgelist));
	FILE *file;

	el->n=0;
	el->e=0;
	file=fopen(input,"r");
	el->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%lu %lu", &(el->edges[el->e].s), &(el->edges[el->e].t))==2) {//Add one edge
		el->n=max3(el->n,el->edges[el->e].s,el->edges[el->e].t);
		el->e++;
		if (el->e==e1) {
			e1+=NLINKS;
			el->edges=realloc(el->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	el->n++;

	el->edges=realloc(el->edges,el->e*sizeof(edge));

	return el;
}

void relabel(edgelist *el){
	Edge i;
	Node source, target, tmp;

	for (i=0;i<el->e;i++) {
		source=el->rank[el->edges[i].s];
		target=el->rank[el->edges[i].t];
		if (source<target){
			tmp=source;
			source=target;
			target=tmp;
		}
		el->edges[i].s=source;
		el->edges[i].t=target;
	}

}

///// CORE ordering /////////////////////

typedef struct {
	Node key;
	Node value;
} keyvalue;

typedef struct {
	Node n_max;	// max number of nodes.
	Node n;	// number of nodes.
	Node *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;


bheap *construct(Node n_max){
	Node i;
	bheap *heap=malloc(sizeof(bheap));

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=malloc(n_max*sizeof(Node));
	for (i=0;i<n_max;i++) heap->pt[i]=-1;
	heap->kv=malloc(n_max*sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap,Node i, Node j) {
	keyvalue kv_tmp=heap->kv[i];
	Node pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

void bubble_up(bheap *heap,Node i) {
	Node j=(i-1)/2;
	while (i>0) {
		if (heap->kv[j].value>heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

void bubble_down(bheap *heap) {
	Node i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->kv[j2].value<heap->kv[j1].value) ) ? j2 : j1 ;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

void insert(bheap *heap,keyvalue kv){
	heap->pt[kv.key]=(heap->n)++;
	heap->kv[heap->n-1]=kv;
	bubble_up(heap,heap->n-1);
}

void update(bheap *heap,Node key){
	Node i=heap->pt[key];
	if (i!=-1){
		((heap->kv[i]).value)--;
		bubble_up(heap,i);
	}
}

keyvalue popmin(bheap *heap){
	keyvalue min=heap->kv[0];
	heap->pt[min.key]=-1;
	heap->kv[0]=heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key]=0;
	bubble_down(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(Node n,Node *v){
	Node i;
	keyvalue kv;
	bheap* heap=construct(n);
	for (i=0;i<n;i++){
		kv.key=i;
		kv.value=v[i];
		insert(heap,kv);
	}
	return heap;
}

void freeheap(bheap *heap){
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//computing degeneracy ordering and core value
void ord_core(edgelist* el){
	Node i, r=0, n=el->n;
	Edge j, e=el->e;
	keyvalue kv;
	bheap *heap;

	Node *d0=calloc(el->n,sizeof(Node));
	Edge *cd0=malloc((el->n+1)*sizeof(Edge));
	Node *adj0=malloc(2*el->e*sizeof(Node));
	for (j=0;j<e;j++) {
		d0[el->edges[j].s]++;
		d0[el->edges[j].t]++;
	}
	cd0[0]=0;
	for (i=1;i<n+1;i++) {
		cd0[i]=cd0[i-1]+d0[i-1];
		d0[i-1]=0;
	}
	for (j=0;j<e;j++) {
		adj0[ cd0[el->edges[j].s] + d0[ el->edges[j].s ]++ ]=el->edges[j].t;
		adj0[ cd0[el->edges[j].t] + d0[ el->edges[j].t ]++ ]=el->edges[j].s;
	}

	heap=mkheap(n,d0);

	el->rank=malloc(n*sizeof(Node));
	for (i=0;i<n;i++){
		kv=popmin(heap);
		el->rank[kv.key]=n-(++r);
		for (j=cd0[kv.key];j<cd0[kv.key+1];j++){
			update(heap,adj0[j]);
		}
	}
	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}

//////////////////////////
//Building the special graph
graph* mkgraph(edgelist *el){
	Node i,max;
	Edge j;
	Node *d;
	graph* g=malloc(sizeof(graph));

	d=calloc(el->n,sizeof(Node));

	for (i=0;i<el->e;i++) {
		d[el->edges[i].s]++;
	}

	g->cd=malloc((el->n+1)*sizeof(Edge));
	g->cd[0]=0;
	max=0;
	for (i=1;i<el->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		max=(max>d[i-1])?max:d[i-1];
		d[i-1]=0;
	}
	printf("core value (max truncated degree) = %lu\n",max);
	fflush(stdout);

	g->adj=malloc(el->e*sizeof(Node));

	for (j=0;j<el->e;j++) {
		g->adj[ g->cd[el->edges[j].s] + d[ el->edges[j].s ]++ ]=el->edges[j].t;
	}

	free(d);
	g->core=max;
	g->n=el->n;

	free(el->rank);
	g->edges=el->edges;
	g->e=el->e;
	free(el);

	return g;
}


subgraph* allocsub(graph *g,Kvalue k){
	Kvalue i;
	Node j;
	subgraph* sg=malloc(sizeof(subgraph));
	sg->n=calloc(k,sizeof(Node));
	sg->d=malloc(k*sizeof(Node*));
	sg->nodes=malloc(k*sizeof(Node*));
	for (i=2;i<k;i++){
		sg->d[i]=malloc(g->core*sizeof(Node));
		sg->nodes[i]=malloc(g->core*sizeof(Node));
	}
	sg->lab=calloc(g->core,sizeof(Kvalue));
	sg->adj=malloc(g->core*g->core*sizeof(Node));
	sg->core=g->core;
	
	sg->new=malloc(g->n*sizeof(Node));
	sg->old=malloc(g->core*sizeof(Node));
	for (j=0;j<g->n;j++){
		sg->new[j]=-1;
	}
	return sg;
}

void mksub(graph* g,edge ed,subgraph* sg,Kvalue k){
	Node i,j,l,x,y;
	Node u=ed.s,v=ed.t;
	Edge iL;

	for (i=0;i<sg->n[k-1];i++){
		sg->lab[i]=0;
	}

	for (iL=g->cd[v];iL<g->cd[v+1];iL++){
		sg->new[g->adj[iL]]=-2;
	}

	j=0;
	for (iL=g->cd[u];iL<g->cd[u+1];iL++){
		x=g->adj[iL];
		if (sg->new[x]==-2){
			sg->new[x]=j;
			sg->old[j]=x;
			sg->lab[j]=k-2;
			sg->nodes[k-2][j]=j;
			sg->d[k-2][j]=0;//new degrees
			j++;
		}
	}

	sg->n[k-2]=j;

	for (i=0;i<sg->n[k-2];i++){//reodering adjacency list and computing new degrees
		x=sg->old[i];
		for (l=g->cd[x];l<g->cd[x+1];l++){
			y=g->adj[l];
			j=sg->new[y];
			if (j<-2){
				sg->adj[sg->core*i+sg->d[k-2][i]++]=j;
			}
		}
	}

	for (iL=g->cd[v];iL<g->cd[v+1];iL++){
		sg->new[g->adj[iL]]=-1;
	}
}

void kclique_thread(Kvalue l, subgraph *sg, Clique *n) {
	Node i,j,k,end,u,v,w;

	if(l==2){
		for(i=0; i<sg->n[2]; i++){//list all edges
			u=sg->nodes[2][i];
			end=u*sg->core+sg->d[2][u];
			for (j=u*sg->core;j<end;j++) {
				(*n)++;//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
			}
		}
		return;
	}

	for(i=0; i<sg->n[l]; i++){
		u=sg->nodes[l][i];
		//printf("%u %u\n",i,u);
		sg->n[l-1]=0;
		end=u*sg->core+sg->d[l][u];
		for (j=u*sg->core;j<end;j++){//relabeling nodes and forming U'.
			v=sg->adj[j];
			if (sg->lab[v]==l){
				sg->lab[v]=l-1;
				sg->nodes[l-1][sg->n[l-1]++]=v;
				sg->d[l-1][v]=0;//new degrees
			}
		}
		for (j=0;j<sg->n[l-1];j++){//reodering adjacency list and computing new degrees
			v=sg->nodes[l-1][j];
			end=sg->core*v+sg->d[l][v];
			for (k=sg->core*v;k<end;k++){
				w=sg->adj[k];
				if (sg->lab[w]==l-1){
					sg->d[l-1][v]++;
				}
				else{
					sg->adj[k--]=sg->adj[--end];
					sg->adj[end]=w;
				}
			}
		}

		kclique_thread(l-1, sg, n);

		for (j=0;j<sg->n[l-1];j++){//restoring labels
			v=sg->nodes[l-1][j];
			sg->lab[v]=l;
		}

	}
}

Clique kclique_main(Kvalue k, graph *g) {
	Edge i;
	Clique n=0;
	subgraph *sg;
	#pragma omp parallel private(sg) reduction(+:n)
	{
		sg=allocsub(g,k);
		#pragma omp for schedule(dynamic,1)// nowait
		for(i=0; i<g->e; i++){
			mksub(g,g->edges[i],sg,k);
			kclique_thread(k-2, sg, &n);
		}
	}
	return n;
}

int main(int argc,char** argv){
	edgelist* el;
	graph* g;
	Kvalue k=atoi(argv[2]);
	Clique n;

	omp_set_num_threads(atoi(argv[1]));

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[3]);
	fflush(stdout);
	
	el=readedgelist(argv[3]);
	printf("Number of nodes = %lu\n",el->n);
	printf("Number of edges = %llu\n",el->e);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Building the graph structure\n");
	fflush(stdout);
	
	ord_core(el);
	relabel(el);
	g=mkgraph(el);

	printf("Number of nodes (degree > 0) = %lu\n",g->n);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Iterate over all cliques\n");
	fflush(stdout);
	
	n=kclique_main(k, g);

	printf("Number of %u-cliques: %llu\n",k,n);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	fflush(stdout);
	t1=t2;

	free_graph(g);

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));
	fflush(stdout);

	return 0;
}
