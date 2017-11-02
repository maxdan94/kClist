/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques. It find a k(1+eps) approximation of the kclique densest. It is highly scallable.

To compile:
"gcc kClistDens.c -O9 -o kClistDens -fopenmp".

To execute:
"./kClistDens p k eps edgelist.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
k is the size of the k-cliques
p is the number of threads
eps is the precision
Will print the number of k-cliques and the density of the found kclique densest.
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <omp.h>


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
	unsigned *rank;//ranking of the nodes according to degeneracy ordering
	//unsigned *map;//oldID newID correspondance NOT USED IN THIS VERSION
} edgelist;

typedef struct {
	unsigned n;
	unsigned e;
	edge *edges;//ading this again here: TO IMPROVE
	unsigned *cd;//cumulative degree: (starts with 0) length=n+1
	unsigned *adj;//truncated list of neighbors
	unsigned core;//core value of the graph
} graph;

typedef struct {
	unsigned *n;//n[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *adj;//truncated list of neighbors
	unsigned char *lab;//lab[i] label of node i
	unsigned **nodes;//sub[l]: nodes in G_l
	unsigned core;
} subgraph;

void free_graph(graph *g){
	free(g->cd);
	free(g->adj);
	free(g);
}

void free_subgraph(subgraph *sg, unsigned char k){
	unsigned char i;
	free(sg->n);
	for (i=1;i<k;i++){
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
inline unsigned int max3(unsigned int a,unsigned int b,unsigned int c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

edgelist* readedgelist(char* input){
	unsigned e1=NLINKS;
	edgelist *el=malloc(sizeof(edgelist));
	FILE *file;

	el->n=0;
	el->e=0;
	file=fopen(input,"r");
	el->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u", &(el->edges[el->e].s), &(el->edges[el->e].t))==2) {//Add one edge
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
	unsigned i, source, target, tmp;
	el->n=0;
		//FILE* file=fopen("debug.txt","w");
	for (i=0;i<el->e;i++) {
		source=el->rank[el->edges[i].s];
		target=el->rank[el->edges[i].t];
		if (source<target){
			tmp=source;
			source=target;
			target=tmp;
		}
		if (source+1>el->n){
			el->n=source+1;
		}
		el->edges[i].s=source;
		el->edges[i].t=target;
		//fprintf(file,"%u %u\n",source,target);
	}
	//fclose(file);

}

///// CORE ordering /////////////////////

typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;


bheap *construct(unsigned n_max){
	unsigned i;
	bheap *heap=malloc(sizeof(bheap));

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=malloc(n_max*sizeof(unsigned));
	for (i=0;i<n_max;i++) heap->pt[i]=-1;
	heap->kv=malloc(n_max*sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap,unsigned i, unsigned j) {
	keyvalue kv_tmp=heap->kv[i];
	unsigned pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

void bubble_up(bheap *heap,unsigned i) {
	unsigned j=(i-1)/2;
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
	unsigned i=0,j1=1,j2=2,j;
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

void update(bheap *heap,unsigned key){
	unsigned i=heap->pt[key];
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
bheap* mkheap(unsigned n,unsigned *v){
	unsigned i;
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
	unsigned i,j,r=0,n=el->n,e=el->e;
	keyvalue kv;
	bheap *heap;

	unsigned *d0=calloc(el->n,sizeof(unsigned));
	unsigned *cd0=malloc((el->n+1)*sizeof(unsigned));
	unsigned *adj0=malloc(2*el->e*sizeof(unsigned));
	for (i=0;i<e;i++) {
		d0[el->edges[i].s]++;
		d0[el->edges[i].t]++;
	}
	cd0[0]=0;
	for (i=1;i<n+1;i++) {
		cd0[i]=cd0[i-1]+d0[i-1];
		d0[i-1]=0;
	}
	for (i=0;i<e;i++) {
		adj0[ cd0[el->edges[i].s] + d0[ el->edges[i].s ]++ ]=el->edges[i].t;
		adj0[ cd0[el->edges[i].t] + d0[ el->edges[i].t ]++ ]=el->edges[i].s;
	}

	heap=mkheap(n,d0);

	el->rank=malloc(n*sizeof(unsigned));
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
	unsigned i,max;
	unsigned *d;
	graph* g=malloc(sizeof(graph));

	d=calloc(el->n,sizeof(unsigned));

	for (i=0;i<el->e;i++) {
		d[el->edges[i].s]++;
	}

	g->cd=malloc((el->n+1)*sizeof(unsigned));
	g->cd[0]=0;
	max=0;
	for (i=1;i<el->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		max=(max>d[i-1])?max:d[i-1];
		d[i-1]=0;
	}
	printf("core value (max truncated degree) = %u\n",max);

	g->adj=malloc(el->e*sizeof(unsigned));

	for (i=0;i<el->e;i++) {
		g->adj[ g->cd[el->edges[i].s] + d[ el->edges[i].s ]++ ]=el->edges[i].t;
	}

	free(d);
	g->core=max;
	g->n=el->n;

	free(el->rank);
	g->edges=el->edges;
	g->e=el->e;
	//free(el);
	//		printf("el2=%u\n",el->e);
	return g;
}


subgraph* allocsub(graph *g,unsigned char k){
	unsigned i;
	subgraph* sg=malloc(sizeof(subgraph));
	sg->n=calloc(k,sizeof(unsigned));
	sg->d=malloc(k*sizeof(unsigned*));
	sg->nodes=malloc(k*sizeof(unsigned*));
	for (i=1;i<k;i++){/////////
		sg->d[i]=malloc(g->core*sizeof(unsigned));
		sg->nodes[i]=malloc(g->core*sizeof(unsigned));
	}
	sg->lab=calloc(g->core,sizeof(unsigned char));
	sg->adj=malloc(g->core*g->core*sizeof(unsigned));
	sg->core=g->core;
	return sg;
}


unsigned *old=NULL,*new=NULL;//to improve
#pragma omp threadprivate(new,old)

void mksub(graph* g,edge ed,subgraph* sg,unsigned char k){
	unsigned i,j,l,x,y;
	unsigned u=ed.s,v=ed.t;

	if (old==NULL){
		new=malloc(g->n*sizeof(unsigned));
		old=malloc(g->core*sizeof(unsigned));
		for (i=0;i<g->n;i++){
			new[i]=-1;
		}
	}

	for (i=0;i<sg->n[k-1];i++){
		sg->lab[i]=0;
	}

	for (i=g->cd[v];i<g->cd[v+1];i++){
		new[g->adj[i]]=-2;
	}

	j=0;
	for (i=g->cd[u];i<g->cd[u+1];i++){
		x=g->adj[i];
		if (new[x]==-2){
			new[x]=j;
			old[j]=x;
			sg->lab[j]=k-2;
			sg->nodes[k-2][j]=j;
			sg->d[k-2][j]=0;//new degrees
			j++;
		}
	}

	sg->n[k-2]=j;

	for (i=0;i<sg->n[k-2];i++){//reodering adjacency list and computing new degrees
		x=old[i];
		for (l=g->cd[x];l<g->cd[x+1];l++){
			y=g->adj[l];
			j=new[y];
			if (j<-2){
				sg->adj[sg->core*i+sg->d[k-2][i]++]=j;
			}
		}
	}

	for (i=g->cd[v];i<g->cd[v+1];i++){
		new[g->adj[i]]=-1;
	}
}


unsigned long long *ckdeg_p,*ckdeg;
unsigned *ck_p;
#pragma omp threadprivate(ckdeg_p,ck_p)

void allocglobal(graph *g,unsigned k){
        #pragma omp parallel
        {
                ck_p=calloc(k,sizeof(unsigned));
                ckdeg_p=calloc(g->n,sizeof(unsigned long long));
        }
        ckdeg=calloc(g->n,sizeof(unsigned long long));
}

void kclique_thread(unsigned char kmax, unsigned char l, subgraph *sg, unsigned long long *n) {
	unsigned i,j,k,end,u,v,w;

	if (kmax==3){//can be improved
		for(i=0; i<sg->n[1]; i++){//list all nodes
			ckdeg_p[old[sg->nodes[1][i]]]++;
			ckdeg_p[ck_p[1]]++;
			ckdeg_p[ck_p[2]]++;
			(*n)++;//listing here!!!
		}
		return;
	}

	if(l==2){
		for(i=0; i<sg->n[2]; i++){//list all edges
			u=sg->nodes[2][i];
			end=u*sg->core+sg->d[2][u];
			for (j=u*sg->core;j<end;j++) {
				ckdeg_p[old[sg->adj[j]]]++;
				ckdeg_p[old[u]]++;
				for (l=2;l<kmax;l++){//ok to use l here :)
					ckdeg_p[ck_p[l]]++;
				}
				(*n)++;//listing here!!!
			}
		}
		return;
	}

	for(i=0; i<sg->n[l]; i++){
		u=sg->nodes[l][i];
		ck_p[l-1]=old[u];
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

		kclique_thread(kmax,l-1, sg, n);

		for (j=0;j<sg->n[l-1];j++){//restoring labels
			v=sg->nodes[l-1][j];
			sg->lab[v]=l;
		}

	}
}

unsigned long long kclique_main(unsigned char k, graph *g) {
	unsigned i;
	unsigned long long n=0;
	subgraph *sg;
	#pragma omp parallel private(sg,i) reduction(+:n)
	{
		sg=allocsub(g,k);

		#pragma omp for schedule(dynamic, 1) nowait
		for(i=0; i<g->e; i++){
			ck_p[k-1]=g->edges[i].s;
			ck_p[k-2]=g->edges[i].t;
			mksub(g,g->edges[i],sg,k);
			kclique_thread(k,k-2, sg, &n);
		}

		free_subgraph(sg,k);

		#pragma omp single
		{
		bzero(ckdeg,g->n*sizeof(unsigned long long));
		}

		#pragma omp barrier //is it necessary???

		#pragma omp critical
		{
			for(i=0; i<g->n; i++){
				ckdeg[i]+=ckdeg_p[i];
			}
			bzero(ckdeg_p,g->n*sizeof(unsigned long long));
		}

	}
	return n;
}

void rmnodes(bool *rm,edgelist* el){
	unsigned long long i,r=0;
	FILE* file=fopen("debug.txt","w");
		for (i=0;i<el->e;i++){
				if (((rm[el->edges[i].s]==1) || (rm[el->edges[i].t]==1)) == 0){
					r++;
									fprintf(file,"%u %u\n",el->edges[i].s,el->edges[i].t);
				}
		}
	//printf("el0=%llu\n",r);
	for (i=0;i<el->e;i++){
//printf("%llu\n",i);
		if ((rm[el->edges[i].s]==1) || (rm[el->edges[i].t]==1)){
			el->edges[i--]=el->edges[--(el->e)];
		}
	}
			//printf("el1=%u\n",el->e);
			fclose(file);
}

int main(int argc,char** argv){
	edgelist* el;
	graph* g;
	unsigned char k=atoi(argv[2]);
	double eps=atof(argv[3]);
	unsigned long long nck;

	omp_set_num_threads(atoi(argv[1]));

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[4]);

	el=readedgelist(argv[4]);
	printf("Number of nodes = %u\n",el->n);
	printf("Number of edges = %u\n",el->e);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Building the graph structure\n");
	ord_core(el);
	relabel(el);
	g=mkgraph(el);

	printf("Number of nodes (degree > 0) = %u\n",g->n);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;


	unsigned i,n_m,e_m;
	double rho,rho_m=0,erho,erho_m;
	unsigned long long nck_m;
	bool *rm=calloc(g->n,sizeof(bool));
	allocglobal(g,k);//allocataing global variables
	nck=kclique_main(k, g);
	do{
		printf("\nCurrent graph:\n");
		printf("Number of %u-cliques: %llu\n",k,nck);
		printf("Number of nodes: %u\n",g->n);
		printf("Number of edges: %u\n",g->e);
		erho=2.*((double)g->e)/((double)(g->n*(g->n-1)));
		printf("edge density: %le\n",erho);
		rho=((double)nck)/((double)(g->n));
		printf("Density: %le\n",rho);
		if (rho>rho_m){
			nck_m=nck;
			rho_m=rho;
			erho_m=erho;
			e_m=g->e;
			n_m=g->n;
		}
		printf("\nDensest found:\n");
		printf("Number of %u-cliques: %llu\n",k,nck_m);
		printf("Number of nodes: %u\n",n_m);
		printf("Number of edges: %u\n",e_m);
		printf("edge density: %le\n",erho_m);
		printf("Density: %le\n",rho_m);
		rho*=k*(1.+eps);
		printf("rho = %e\n",rho);
		//unsigned long long r=0,r2=0;
		//printf("n,m,n,m= %u %u %u %u\n",g->n,g->e,el->n,el->e);
		for (i=0;i<g->n;i++){
			//if (rm[i]==1){
			//	printf("aaaaaaaaaaa\n");
			//}
			//r+=ckdeg[i];
			if (((double)ckdeg[i])<rho){
				rm[i]=1;
			}
			//else{
			//	printf("not remove: %u, %llu\n",i,ckdeg[i]);
			//					r2++;
			//}
		}
		//printf("sumckdeg, nnodes= %llu, %llu\n",r,r2);
		rmnodes(rm,el);
		bzero(rm,g->n*sizeof(bool));
		ord_core(el);
		relabel(el);
			//printf("n,m= %u %u\n",el->n,el->e);
		free_graph(g);
		g=mkgraph(el);
		nck=kclique_main(k, g);
		//printf("nck=%llu\n",nck);
	}while (nck>0);


	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	printf("\n%u %u %u %llu %e %e\n",k,n_m,e_m,nck_m,erho_m,rho_m);

	return 0;
}
