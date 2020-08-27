#include <stdio.h>
#include <stdlib.h>
#define MAX_V 64
#define BOARD_SIZE 8
/********* TAD GRAPH *******/
typedef struct t_graph {
    int V;
    int E;
    int **adj;
    int *degree;
    int max_degree;
}Graph;
typedef struct t_graph *TGraph;

typedef struct t_edge{
    int u;
    int v;
}Edge;

int* Array_init(int size){
    int *array = calloc(size, sizeof(int));
    
    return array;
}

int** Matrix_init(int row_num, int col_num, int value){
    int i, j;

    int **matrix = malloc(row_num*sizeof(int*));

    for(i = 0; i < row_num; i++){
        matrix[i] = malloc(col_num*sizeof(int));
        
        for(j = 0; j < col_num; j++)
            matrix[i][j] = value;
    }

    return matrix;    
}

TGraph Graph_init(int V, int max_degree, int value){
    TGraph graph = (TGraph)malloc(sizeof(Graph));

    if(graph != NULL){
        
        graph->V = V;
        graph->E = 0;
        graph->max_degree = max_degree;
        graph->degree = Array_init(V);
        graph->adj = Matrix_init(MAX_V, max_degree, 0);
    }

    return graph;
}

void Graph_destroy(TGraph graph){
    if(graph != NULL){
        for (int i = 0; i < graph->V; i++)
            free(graph->adj[i]);
        free(graph->adj);
        free(graph->degree);
        free(graph);
    }
}

int Graph_insertE(TGraph graph, Edge e){
    if(graph==NULL) return 0;
    if(e.u < 0 || e.u >= graph->V) return 0;
    if(e.v < 0 || e.v >=graph->V) return 0;

    graph->adj[e.u][graph->degree[e.u]] = e.v;
    graph->degree[e.u]++;
    graph->adj[e.v][graph->degree[e.v]] = e.u;
    graph->degree[e.v]++;

    graph->E++;
    
    return 1;
}

int Graph_removeE(TGraph graph, Edge e){
    if(graph==NULL) return 0;
    if(e.u < 0 || e.u >= graph->V) return 0;
    if(e.v < 0 || e.v >=graph->V) return 0;
    
    int i=0;
    while (i<graph->degree[e.u] && graph->adj[e.u][i]!= e.v) 
        i++;
    
    if(i==graph->degree[e.u]) return 0;

    graph->degree[e.u]--;
    graph->E--;
    graph->adj[e.u][i] = graph->adj[e.u][graph->degree[e.u]];

    i=0;
    while (i<graph->degree[e.v] && graph->adj[e.v][i]!= e.u) 
        i++;
    if(i==graph->degree[e.v]) return 0;
    graph->degree[e.u]--;
    graph->adj[e.v][i] = graph->adj[e.v][graph->degree[e.v]];

    return 1;
}

void Graph_print(TGraph graph){
    int i, j;
    for(i=0; i<graph->V; i++){
        printf("[%d]->", i);
        for(j=0; j<graph->degree[i]; j++){
             printf("%d ", graph->adj[i][j]);    
        }
        printf("\n");
    }
}

int Graph_findEdge(TGraph graph, Edge e){
    int i;

    for(i=0; i<graph->degree[e.u]; i++){
        if(graph->adj[e.u][i]==e.v) return 1;
    }
    return 0;
}

void Edge_print(Edge e){
    printf("edge(%d, %d)\n", e.u, e.v);
}

void Array_print(int *v, int size){
    int i;
    for ( i = 0; i < size; i++)
    {
        printf("%d ", v[i]);
    }
    printf("\n");
    
}
/********Knight Path***********/
void Verify_n_insert(TGraph graph, Edge e){
    int value, col_dist, row_dist;
    //calculates distance between svertices
    col_dist = abs(e.u%BOARD_SIZE - e.v%BOARD_SIZE);
    row_dist = abs((int)e.u/BOARD_SIZE - (int)e.v/BOARD_SIZE);
    value = col_dist + row_dist;
    //Edge_print(e);
    //printf("edge(%d, %d) ------ value = %d\n", e.u, e.v, value);
    if(value == 3 && col_dist != 3 && row_dist != 3 && !Graph_findEdge(graph, e)) Graph_insertE(graph, e);
}

void KnightMove_insert(TGraph graph){
    int i,j,k;
    Edge e;
    for(i=0; i<graph->V; i++){
        e.u = i;
        for(j=0; j<graph->V; j++){
            e.v = j;
            Verify_n_insert(graph, e);
        }
    }
    //return graph;
}

int findMinDist(int *dist, int *visited, int V){
    int i, min = -1, first = 1;

    for(i=0; i< V; i++){
        if(dist[i]>=0 && visited[i]==0){
            if(first){
                min = i;
                first = 0;
            }
            else{
                if(dist[min]>dist[i])
                    min = i;
            }
        }
    }
    return min;
}

void Graph_minPath(TGraph graph, Edge e, int *path, int *dist){
    int i, count, ind, *visited, aux;

    count = graph->V;
    visited = (int*) malloc(graph->V*sizeof(int));
    for ( i = 0; i < graph->V; i++)
    {
        path[i]=-1;
        dist[i]=-1;
        visited[i]=0;
    }
    dist[e.u] = 0;

    while (count >0)
    {
        aux = findMinDist(dist, visited, graph->V);
        if(aux == -1) break;

        visited[aux] = 1;

        count --;

        for ( i = 0; i < graph->degree[aux]; i++){
            ind = graph->adj[aux][i];

            if(dist[ind]<0){
                dist[ind] = dist[aux]+1;
                path[ind] = aux;
            }
            else{
                if(dist[ind] > dist[aux]+1){
                    dist[ind] = dist[aux];
                    path[ind] = aux;
                }
            }
        }
    }

    free(visited);
}

int main()
{
    //input example "c1 b5" <=> "initial position  final position" of knight
    char cci, ccf; //char columns identificator a-h
    int ci, cf, ri, rf; //columns and rows of initial position and final position
    Edge e; //this edge will store the initial position and final position

    scanf("%c%d %c%d", &cci, &ri, &ccf, &rf);
    
    //converts char column in int column identifier
    ci = cci - 97;
    cf = ccf - 97;
    
    //calculates initial and final vertices
    e.u = (ri-1)*8+ci;
    e.v = (rf-1)*8+cf;

    //Edge_print(e);
    TGraph graph;
    graph = Graph_init(MAX_V, 8, 0);
    KnightMove_insert(graph);

    int dist[MAX_V], path[MAX_V];

    //Graph_print(graph);
    
    Graph_minPath(graph, e, path, dist);
    //Array_print(path, MAX_V);
    printf("%d\n", dist[e.v]);

    Graph_destroy(graph);
    return 0;
}
