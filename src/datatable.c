/*
MIT License

Copyright (c) 2019 Gleb Goussarov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <string.h>
#include "vecops.h"
#include "datatable.h"
#include "datamap.h"
#include "stats.h"

struct datatable_t {
    double** values;
    size_t w;
    size_t h;
    DM64_t* rowids;
    DM64_t* colids;
    char** rownames;
    char** colnames;
    double fill;
};

datatable_t* datatable_alloc(size_t ncol, size_t nrow, double autofill) {
    datatable_t* result;
    size_t i,j;
    result = (datatable_t*) calloc(1, sizeof(datatable_t));
    result->w = ncol;
    result->h = nrow;
    result->values = (double**)calloc(nrow, sizeof(double*));
    for (i = 0;i < nrow;i++) {
        result->values[i] = (double*)malloc(ncol*sizeof(double));
        for (j = 0;j < ncol;j++) {
            result->values[i][j] = autofill;
        }
    }
    result->rowids = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
    result->colids = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
    result->rownames = (char**)calloc(nrow, sizeof(char*));
    result->colnames = (char**)calloc(ncol, sizeof(char*));
    result->fill = autofill;
    return result;
}
void datatable_free(datatable_t* target) {
    size_t i;
    for (i = 0;i < target->h;i++) {
        free(target->values[i]);
        target->values[i] = NULL;
        if (target->rownames[i])free(target->rownames[i]);
        target->rownames[i] = NULL;
    }
    free(target->values);
    free(target->rownames);
    for (i = 0;i < target->w;i++) {
        if (target->colnames[i])free(target->colnames[i]);
        target->colnames[i] = NULL;
    }
    free_DM64(target->colids);
    free_DM64(target->rowids);
    target->colids = NULL;
    target->rowids = NULL;
}

size_t datatable_ncols(datatable_t* target) {
    return target->w;
}
size_t datatable_nrows(datatable_t* target){
    return target->h;
}

void _datatable_addcolsupto(datatable_t* target, size_t col) {
    size_t i;
    size_t j;
    if (col >= target->w) {
        for (i = 0;i < target->h;i++) {
            target->values[i] = (double*)realloc(target->values[i], sizeof(double)*(col + 1));
            for (j = target->w; j <= col;j++) {
                target->values[i][j] = target->fill;
            }
        }
        target->colnames = realloc(target->colnames, sizeof(char*)*(col + 1));
        for (j = target->w; j <= col;j++) {
            target->colnames[j] = NULL;
        }
        target->w = col + 1;
    }
}
void _datatable_addrowsupto(datatable_t* target, size_t row) {
    size_t i;
    size_t j;
    if(row >= target->h) {
        target->values = (double**)realloc(target->values, sizeof(double*)*(row + 1));
        target->rownames = (char**)realloc(target->rownames, sizeof(char*)*(row + 1));
        for (i = target->h;i <= row;i++) {
            target->values[i] = (double*)malloc(sizeof(double)*target->w);
            target->rownames[i] = NULL;
            for (j = 0; j < target->w;j++) {
                target->values[i][j] = target->fill;
            }
        }
        target->h = row + 1;
    }
}

errcode_t datatable_set(datatable_t* target, size_t col, size_t row, double value) {
    _datatable_addcolsupto(target, col);
    _datatable_addrowsupto(target, row);
    target->values[row][col] = value;
    return 0;
}
double datatable_get(datatable_t* target, size_t col, size_t row, int* nullflag){
    if (col >= target->w || row >= target->h) {
        *nullflag = 1;
        return target->fill;
    }
    return target->values[row][col];
}

void datatable_addtables(datatable_t* target, datatable_t** toadd, size_t numtables) {
    size_t i,j,k;
    double newvalue;
    size_t targcol;
    size_t targrow;
    char* name;
    int nf;
    for (i = 0;i < numtables;i++) {
        for (j = 0;j < datatable_nrows(toadd[i]);j++) {
            name = datatable_getrowname(toadd[i], j);
            if (name)
                targrow = datatable_getrowid(target, name);
            else
                targrow = j;
            for (k = 0;k < datatable_ncols(toadd[i]);k++) {
                name = datatable_getcolname(toadd[i], k);
                if (name)
                    targcol = datatable_getcolid(target, name);
                else
                    targcol = k;
                newvalue = datatable_get(target, targcol, targrow, &nf);
                newvalue += datatable_get(toadd[i], k, j, &nf);
                datatable_set(target,targcol, targrow, newvalue);
            }
        }
    }
}

void datatable_setcolname(datatable_t* target, size_t col, const char* name) {
    size_t namelen;
    _datatable_addcolsupto(target, col);
    if (target->colnames[col]) free(target->colnames[col]);
    namelen = strlen(name) + 1;
    target->colnames[col] = malloc(namelen);
    memcpy(target->colnames[col], name, namelen);
    DM64_append(target->colids, name, (int)namelen - 1, (int64_t)col);
}
size_t datatable_getcolid(datatable_t* target, const char* name){
    int nullflag;
    int64_t result;
    result = DM64_get(target->colids, name, (int) strlen(name), &nullflag);
    if (nullflag)
        result = target->w;
    return result;
}
char* datatable_getcolname(datatable_t* target, size_t col) {
    if (col >= target->w)return NULL;
    return target->colnames[col];
}

void datatable_setrowname(datatable_t* target, size_t row, const char* name) {
    size_t namelen;
    _datatable_addrowsupto(target, row);
    if (target->rownames[row]) free(target->rownames[row]);
    namelen = strlen(name) + 1;
    target->rownames[row] = malloc(namelen);
    memcpy(target->rownames[row], name, namelen);
    DM64_append(target->rowids, name, (int)namelen - 1, (int64_t)row);
}
ssize_t datatable_getrowid(datatable_t* target, const char* name) {
    int nullflag;
    int64_t result;
    result = DM64_get(target->rowids, name, (int) strlen(name), &nullflag);
    if (nullflag)
        result = target->h;
    return result;
}
char* datatable_getrowname(datatable_t* target, size_t row) {
    if (row >= target->h)return NULL;
    return target->rownames[row];
}

size_t _find_medoid(double** data, size_t begin, size_t end, int64_t* curorder) {
    double sumofdists;
    size_t result;
    size_t i, j;
    double minsod;
    minsod = -1;
    /**/
    for (i = begin;i < end;i++) {
        sumofdists = 0;
        for (j = begin;j < end;j++) {
            sumofdists += data[curorder[i]][curorder[j]];
        }
        if (minsod < 0) {
            minsod = sumofdists;
            result = i;
        }
        else if (minsod >= sumofdists) {
            minsod = sumofdists;
            result = i;
        }
    }
    return result;
}
size_t _distance_matrix_progressivesort(double** data, size_t begin, size_t end, double defaultcut, int64_t** curorder) {
    size_t i, i0, nclusters, j0, tmppivot;
    int64_t* order;
    int64_t tmpi64;
    double* currow;
    double delta;
    double cutat;
    size_t size;
    
    size = end;
    order = *curorder;
    if (!order) {
        order = (int64_t*)malloc(sizeof(int64_t)*size);
        for (i = 0;i < size;i++) order[i] = (int64_t)(i);
    }

    currow = (double*)malloc(sizeof(double)*size);
    
    i0 = begin;
    nclusters = 0;
    j0 = begin;
    cutat = defaultcut;
    while (i0 < end-1) {
        for (i = i0;i < end;i++) {
            currow[i] = data[order[i0]][order[i]];
        }
        vec_sort_by(order + i0, sizeof(int64_t), currow + i0, end - i0);
        delta = data[order[i0]][order[i0 + 1]];
        if (delta > cutat) {
            tmppivot = _find_medoid(data, j0, i0+1, order);
            if (tmppivot != i0) {
                tmpi64 = order[tmppivot];
                order[tmppivot] = order[i0];
                order[i0] = tmpi64;
                for (i = j0;i < i0;i++) {
                    if (data[order[i]][order[i0]] > cutat)cutat = data[order[i]][order[i0]];
                }
            }
            else {
                cutat = defaultcut;
                nclusters++;
                i0++;
                j0 = i0;
            }
        }
        else
            i0++;
    }
    nclusters++;
    free(currow);
    
    *curorder = order;
    return  nclusters;
}
void closest_points(double** distance_matrix, size_t begin, size_t end, int64_t* curorder, size_t** p_result, double** p_dist) {
    size_t* result;
    size_t size;
    size_t i, j, x, y;
    double mindist;
    double curdist;
    double* dist;

    size = end - begin;

    result = *p_result;
    if (!result)
        result = (size_t*)malloc(sizeof(size_t)*size);
    dist = *p_dist;
    if (!dist)
        dist = (double*)malloc(sizeof(double)*size);
    for (i = 0;i < size;i++) {
        y = curorder[i + begin];
        mindist = -1;
        for (j = 0;j < size;j++) {
            x = curorder[j + begin];
            curdist = distance_matrix[y][x];
            if (x != y && (mindist < 0 || mindist > curdist)) {
                mindist = curdist;
                result[i] = j;
            }
        }
        dist[i] = mindist;
    }
    *p_dist = dist;
    *p_result = result;
}

#define DENDRO_INTERLACED_CHILDREN 1
#define DENDRO_ALL_NODES_NAMED 2
void _print_dendrogram_as_newick_tree(int64_t* parent, int64_t* child1, int64_t* child2, double* branch_lengths, double* node_values, char** node_names, size_t size, PF_t* f, uint32_t flags) {
    int64_t i, k;
    size_t child_stride;
    int64_t isize;
    size_t node_offset;
    int64_t root_id;
    int all_nodes_named;
    double branch_length;
    double node_value;
    if (!f)return;
    isize = (int64_t)size;
    if (flags&DENDRO_INTERLACED_CHILDREN) child_stride = 2;
    else child_stride = 1;
    if (flags&DENDRO_ALL_NODES_NAMED) all_nodes_named = 1;
    else all_nodes_named = 0;
    root_id = size * 2 - 2;
    k = root_id;
    /* move down to the "leftmost" leaf */
    while (k >= isize) {
        PFputc('(', f); /* add '(' whenever we go down a level*/
        node_offset = (size_t)(k - isize) * child_stride;
        k = child1[node_offset];
    }
    for (i = 0; i < isize;i++) {
        /*print the name of the current node */
        if (branch_lengths) branch_length = branch_lengths[k];
        else branch_length = 1;
        /* node values are not printed for leaves */
        if (node_names)PFprintf(f, "%s:%f", node_names[k], branch_length);
        else PFprintf(f, "N" _LLD_ ":%f", k, branch_length);
        /*move up until the current node is no longer the last child of its parent*/
        node_offset = (size_t)(parent[k] - size)*child_stride;
        while (k != root_id && k == child2[node_offset]) {
            k = parent[k];
            node_offset = (size_t)(parent[k] - size)*child_stride;
            /* add ')' and generate the node label whenever we go up a level*/
            if (branch_lengths) branch_length = branch_lengths[k];
            else branch_length = 1;
            if (node_values) node_value = node_values[k];
            else node_value = 0.0;

            if (all_nodes_named)
                PFprintf(f, ")%s-%f:%f", node_names[k], node_value, branch_length);
            else
                PFprintf(f, ")C" _LLD_ "-%f:%f", k, node_value, branch_length);
        }
        if (k != root_id) {
            PFputc(',', f); /* add ',' whenever we stay on the same level*/
            k = child2[node_offset];
            while (k >= isize) {
                PFputc('(', f); /* add '(' whenever we go down a level*/
                node_offset = (size_t)(k - isize) * child_stride;
                k = child1[node_offset];
            }
        }
    }
    PFputc(';', f);
}
int64_t _find_side_child(int64_t* child1, int64_t* child2, char* clusterdir, int64_t clustid, int64_t datasize, char side) {
    char curdir;
    curdir = clusterdir[clustid];
    while (clustid >= datasize) {
        if (curdir == side)clustid = child1[clustid - datasize];
        else clustid = child2[clustid - datasize];
        if (clusterdir[clustid]) curdir = (!(curdir));
    }
    return clustid;
}
double _breakpoint_mse(double* function, size_t nelts, size_t breakpoint) {
    double mse;
    double angle;
    double delta;
    double ymid;
    size_t i;
    mse = 0;

    ymid = function[breakpoint];
    angle = ymid / ((double)breakpoint);
    for (i = 0;i < breakpoint;i++) {
        delta = function[i] - angle*(double)i;
        mse += delta*delta;
    }
    angle = (function[nelts - 1] - ymid) / ((double)(nelts - breakpoint));
    for (i = breakpoint;i < nelts;i++) {
        delta = function[i] - (ymid + angle*(double)(i - breakpoint));
        mse += delta*delta;
    }
    return mse;
}
size_t _find_breakpoint(double* function, size_t nelts){
    /* complexity of this function: N*log(N) */
    double oldmse, newmseleft, newmseright;
    size_t left, right;
    size_t breakpoint;

    breakpoint = nelts / 2;
    left = 0;
    right = nelts;

    oldmse = _breakpoint_mse(function, nelts, breakpoint);
    newmseleft = _breakpoint_mse(function, nelts, (left + breakpoint) / 2);
    newmseright = _breakpoint_mse(function, nelts, (right + breakpoint) / 2);

    while (right - left > 4) {
        if (newmseleft < newmseright) {
            if (newmseleft < oldmse) {
                right = breakpoint;
                breakpoint = (left + breakpoint) / 2;
                oldmse = newmseleft;
            }
            else {
                left = (left + breakpoint) / 2;
                right = (right + breakpoint) / 2;
            }
        }
        else {
            if (newmseright < oldmse) {
                left = breakpoint;
                breakpoint = (right + breakpoint) / 2;
                oldmse = newmseright;
            }
            else {
                left = (left + breakpoint) / 2;
                right = (right + breakpoint) / 2;
            }
        }
        newmseleft = _breakpoint_mse(function, nelts, (left + breakpoint) / 2);
        newmseright = _breakpoint_mse(function, nelts, (right + breakpoint) / 2);
    }

    return breakpoint;
}
size_t _distance_matrix_partsort(double** data, size_t begin, size_t end, double sort_at, int64_t** curorder, double* avgtheshold) {
    size_t i, i0, firsti, nclusters, pivot;
    int64_t* order;
    double* currow;
    double curclustlimit;
    double outnext_avg;
    double rel_to_ext;
    double factor;
    size_t size;

    size = end;
    order = *curorder;
    if (!order) {
        order = (int64_t*)malloc(sizeof(int64_t)*size);
        for (i = 0;i < size;i++) order[i] = (int64_t)(i);
    }

    currow = (double*)malloc(sizeof(double)*size);

    i0 = begin;
    nclusters = 0;
    firsti = i0;
    curclustlimit = sort_at;
    outnext_avg = 0;
    while (i0 < end) {
        for (i = firsti;i < end;i++) {
            currow[i] = data[order[i0]][order[i]];
        }
        vec_sort_by(order + i0, sizeof(int64_t), currow + i0, end - i0);

        i = i0;
        while (i<end && currow[i] < curclustlimit)i++;
        if (i - i0 == 1) {
            nclusters++;
            i0 = i;

            /* sort everything according to the medoid, then sort the inside of the cluster according to a near-outlier */
            pivot = _find_medoid(data, firsti, i, order);
            for (i = firsti;i < end;i++) {
                currow[i] = data[order[pivot]][order[i]];
            }
            vec_sort_by(order + firsti, sizeof(int64_t), currow + firsti, end - firsti);
            pivot = i0 - 1;
            if (i0 < end) {
                /* prepare reverse sort according to the point which is closest to the next cluster */
                for (i = firsti;i < i0;i++) {
                    currow[i] = data[order[i0]][order[i]] ;
                }
            }
            else {
                /* prepare reverse sort according to the most unusual point within the cluster */
                for (i = firsti;i < i0;i++) {
                    currow[i] = -data[order[pivot]][order[i]];
                }
            }
            vec_sort_by(order + firsti, sizeof(int64_t), currow + firsti, i0 - firsti);

            rel_to_ext = 0;
            factor = 0;
            if (i0 < end) {
                for (i = firsti;i < i0;i++) {
                    rel_to_ext += data[order[i0]][order[i]];
                }
                factor += 1.0;
            }
            if (firsti > 0) {
                for (i = firsti;i < i0;i++) {
                    rel_to_ext += data[order[firsti-1]][order[i]];
                }
                factor += 1.0;
            }
            if (factor > 0) outnext_avg += rel_to_ext / (((double)(i0 - firsti))*factor);
            firsti = i0;
            curclustlimit = sort_at;
        }
        else {
            i0 = i - 1;
            pivot = _find_medoid(data, firsti, i, order);
            for (i = firsti;i <= i0;i++) {
                currow[i] = data[order[pivot]][order[i]];
            }
            curclustlimit = vec_max(currow + firsti, i0 + 1 - firsti);
            if (curclustlimit < sort_at) curclustlimit = sort_at;
        }
    }
    nclusters++;
    outnext_avg = outnext_avg / ((double)nclusters);
    *avgtheshold = outnext_avg;
    free(currow);

    *curorder = order;
    return  nclusters;
}
void _distance_matrix_slcsort_makedendroandsort(size_t size,double**data, double* L, int64_t* P, int64_t* O, PF_t* output) {
    int64_t* clustalias;
    int64_t* parent;
    int64_t* child1;
    int64_t* child2;
    int64_t* nO;
    char* clusterdir;
    double* branch_length;
    int64_t next, isize, maxid;
    size_t j,j2;
    int64_t cparent, curid, left, right, tmpi64, compareto;
    int64_t lchildleft, lchildright, rchildleft, rchildright;
    double D[4];

    isize = (int64_t)size;

    /* create the dendrogram */
    parent = (int64_t*)malloc(sizeof(int64_t)*size * 2);
    child1 = (int64_t*)malloc(sizeof(int64_t)*size);
    child2 = (int64_t*)malloc(sizeof(int64_t)*size);
    clusterdir = (char*)malloc(sizeof(char)*size * 2);
    clustalias = (int64_t*)malloc(sizeof(int64_t)*size * 2);
    branch_length = (double*)malloc(sizeof(double)*size * 2);
    for (j = 0;j < size;j++) clustalias[j] = (int64_t)j;
    for (j = 0;j < size * 2;j++) parent[j] = -1;
    memset(clusterdir, 0, size * 2);
    next = (int64_t)size;
    for (j = 0;j < size - 1;j++) {
        j2 = j * 2;
        left = clustalias[P[j2 + 1]];
        while (parent[left] != -1) left = parent[left];
        right = clustalias[P[j2]];
        while (parent[right] != -1) right = parent[right];

        lchildleft = _find_side_child(child1, child2, clusterdir, left, isize, 0);
        lchildright = _find_side_child(child1, child2, clusterdir, left, isize, 1);
        rchildleft = _find_side_child(child1, child2, clusterdir, right, isize, 0);
        rchildright = _find_side_child(child1, child2, clusterdir, right, isize, 1);
        D[0] = data[O[lchildleft]][O[rchildleft]];
        D[1] = data[O[lchildright]][O[rchildleft]];
        D[2] = data[O[lchildleft]][O[rchildright]];
        D[3] = data[O[lchildright]][O[rchildright]];
        if (D[1] > D[0]) {
            if (D[0] <= D[2]) {
                if (D[0] <= D[3]) { /* D0 is the smallest distance */
                    clusterdir[left] = (!(clusterdir[left]));
                }
                else { /* D3 is the smallest distance */
                    clusterdir[right] = (!(clusterdir[right]));
                }
            }
            else {
                if (D[2] <= D[3]) { /* D2 is the smallest distance */
                    tmpi64 = left;
                    left = right;
                    right = tmpi64;
                }
                else { /* D3 is the smallest distance */
                    clusterdir[right] = (!(clusterdir[right]));
                }
            }
        }
        else {
            if (D[1] <= D[2]) {
                /* If D1 is the smallest distance, then there is no need to change anything */
                if (D[1] > D[3]) { /* D3 is the smallest distance */
                    clusterdir[right] = (!(clusterdir[right]));
                }
            }
            else {
                if (D[2] <= D[3]) { /* D2 is the smallest distance */
                    tmpi64 = left;
                    left = right;
                    right = tmpi64;
                }
                else { /* D3 is the smallest distance */
                    clusterdir[right] = (!(clusterdir[right]));
                }
            }
        }

        parent[left] = next;
        parent[right] = next;
        child1[next - isize] = left;
        child2[next - isize] = right;
        clustalias[P[j2]] = next;
        clustalias[P[j2 + 1]] = next;
        branch_length[left] = L[j];
        branch_length[right] = L[j];
        branch_length[next] = 0.0;
        next++;
    }
    maxid = next;
    free(clustalias);
    /* reorganize the dendrogram */
    curid = child1[maxid - isize - 1]; /* note that the direction of the last node to be created will always be 0 */
    while (curid >= 0 && curid != maxid - 1) {
        /* move down */
        while (curid >= isize) {
            if (clusterdir[curid] != 0) {
                clusterdir[child1[curid - isize]] = (!(clusterdir[child1[curid - isize]]));
                clusterdir[child2[curid - isize]] = (!(clusterdir[child2[curid - isize]]));
                tmpi64 = child1[curid - isize];
                child1[curid - isize] = child2[curid - isize];
                child2[curid - isize] = tmpi64;
            }
            curid = child1[curid - isize];
        }
        /* move up */
        while (curid != maxid - 1 && curid == child2[parent[curid] - isize]) {
            curid = parent[curid];
        }
        /* move down once if not at the last iteration */
        if (curid >= 0 && curid != maxid - 1) {
            curid = child2[parent[curid] - isize];
        }
    }
    /* print the dendrogram if output is specified */
    _print_dendrogram_as_newick_tree(parent, child1, child2, branch_length, branch_length, NULL, size, output, 0);
    /* sort the order to fit the dendrogram */
    nO = (int64_t*)malloc(sizeof(int64_t)*size);
    curid = maxid - 1;
    while (curid >= isize) {
        curid = child1[curid - isize];
    }
    for (next = 0;next < isize;next++) {
        nO[next] = curid;
        cparent = parent[curid];
        compareto = child2[cparent - isize];
        while (compareto == curid && cparent < maxid) {
            curid = cparent;
            if (curid == maxid - 1) {
                cparent++;
            }
            else cparent = parent[curid];
            compareto = child2[cparent - isize];
        }
        if (cparent < maxid) {
            curid = child2[cparent - isize];
            while (curid >= isize) {
                curid = child1[curid - isize];
            }
        }
    }
    free(parent);
    free(child1);
    free(child2);
    free(clusterdir);
    free(branch_length);
    /* replace local indices with the actual indices */
    for (j = 0;j < size;j++) nO[j] = O[nO[j]];
    for (j = 0;j < size;j++) O[j] = nO[j];
    free(nO);
}
void _distance_matrix_SLINKclust(double** data, int64_t* O, int64_t* P, double* L, size_t size) {
    /* based on the SLINK algorithm described in R. Sibson, "SLINK: And optimally efficient algorithm for the single-link cluster method", The Computer Journal, 1972 */
    /* input:
         data[size][size] - contains dissimilarities between points
         O[size] (order) - defines which data should be looked at
         P[2*size] (Pi) - output for the merging process, contains contains the clust pairs to be merged
         L[size] (Lambda) - output for the merge distances, contains the distance between the clusters to be merged at each step (P)
         size - the number of data points
    */
    double* M;
    size_t i,j,j2;
    M = (double*)malloc(sizeof(double)*size);

    P[0] = 0;
    L[0] = -1.;
    /* clustering using SLINK */
    for (i = 1;i < size;i++) {
        P[i * 2] = i;
        L[i] = -1.;
        for (j = 0;j < i;j++) {
            M[j] = data[O[i]][O[j]];
        }
        for (j = 0;j < i;j++) {
            j2 = j * 2;
            if (L[j] >= M[j] || L[j] < 0) {
                if (M[P[j2]] > L[j] && L[j] > 0) M[P[j2]] = L[j];
                L[j] = M[j];
                P[j2] = i;
            }
            else if (M[P[j2]] > M[j])
                M[P[j2]] = M[j];
        }
        for (j = 0;j < i;j++) {
            if (L[j] >= L[P[j2]] || L[j] < 0) P[j2] = i;
        }
    }
    for (j = 0;j < size;j++) P[j * 2 + 1] = (int64_t)j;
    free(M);
}
size_t _distance_matrix_slcsort(double** data, size_t size, int64_t** order) {
    int64_t* P;
    double* L;
    int64_t* O;
    size_t i;
    size_t break_at;
    size_t nclust;
    double break_value;
    
    if (!order)return 0;
    O = *order;
    if (!O) {
        O = malloc(sizeof(int64_t)*size);
        for (i = 0;i < size;i++) O[i] = (int64_t)i;
    }
    
    P = (int64_t*)malloc(sizeof(int64_t)*size*2);
    L = (double*)malloc(sizeof(double)*size);
    _distance_matrix_SLINKclust(data, O, P, L, size);
    /* sort the cluster joining order using dissimilarity */
    L[size - 1] = vec_max(L, size)+1;
    vec_sort_by(P, sizeof(int64_t)*2, L, size);
    break_at = _find_breakpoint(L, size);
    break_value = L[break_at];
    /*_distance_matrix_slcsort_makedendroandsort(size, data, L, P, O);*/
    free(P);
    free(L);
    *order = O;
    nclust = _distance_matrix_partsort(data,0,size,break_value/2,order,&break_value);
    return nclust;
}
size_t _distance_matrix_slclust(double** data, size_t size, int64_t** order, PF_t* output_f) {
    int64_t* P;
    double* L;
    int64_t* O;
    size_t i;
    size_t break_at;
    size_t nclust;
    double break_value;

    if (!order)return 0;
    O = *order;
    if (!O) {
        O = malloc(sizeof(int64_t)*size);
        for (i = 0;i < size;i++) O[i] = (int64_t)i;
    }

    P = (int64_t*)malloc(sizeof(int64_t)*size * 2);
    L = (double*)malloc(sizeof(double)*size);
    _distance_matrix_SLINKclust(data, O, P, L, size);
    /* sort the cluster joining order using dissimilarity */
    L[size - 1] = vec_max(L, size) + 1;
    vec_sort_by(P, sizeof(int64_t) * 2, L, size);
    break_at = _find_breakpoint(L, size);
    break_value = L[break_at];
    if (break_value > 0) {
        for (i = 0;i < size;i++) {
            L[i] = L[i] / break_value;
        }
    }
    nclust = size - break_at;
    _distance_matrix_slcsort_makedendroandsort(size, data, L, P, O,output_f);
    free(P);
    free(L);
    *order = O;
    return nclust;
}
size_t _distance_matrix_roughclust(double** data, size_t begin, size_t end, double defaultcut, double minjump, int64_t** curorder) {
    size_t i, i0, nclusters, njumps, j;
    int64_t* order;
    int64_t* localorder;
    double* currow;
    double* crdiff;
    double* intradiff;
    double maxdiff, avgdiff, seconddiff;
    double curminjump, oldpivotvalue;
    size_t size;
    size_t pivot, tmppivot, oldpivot;
    size_t skip, maxdiffloc;
    int64_t tmpi64;
    int steptype;

    size = end;
    order = *curorder;
    if (!order) {
        order = (int64_t*)malloc(sizeof(int64_t)*size);
        for (i = 0;i < size;i++) order[i] = (int64_t)(i);
    }
    
    intradiff = (double*)malloc(sizeof(double)*size);
    currow = (double*)malloc(sizeof(double)*size);
    crdiff = (double*)malloc(sizeof(double)*(size-1));
    localorder = (int64_t*)malloc(sizeof(int64_t)*size);

    i0 = begin;
    nclusters = 0;
    skip = 1;
    steptype = 1;
    pivot = i0+1;
    while (i0 < end - skip) {
        for (i = i0;i < end;i++) {
            currow[i] = data[order[i0]][order[i]];
        }
        oldpivotvalue = currow[pivot-1];
        oldpivot = pivot;
        vec_sort_by(order + i0, sizeof(int64_t), currow + i0, end - i0);
        maxdiff = 0;
        maxdiffloc = end;
        seconddiff = 0;
        for (i = i0;i < end - skip;i++) {
            crdiff[i] = currow[i + skip] - currow[i];
            if (crdiff[i] > maxdiff){
                seconddiff = maxdiff;
                maxdiff = crdiff[i];
                maxdiffloc = i;
            }
        }
        avgdiff = 0;
        njumps = 0;
        if (end - i0 > 4) {
            for (i = i0;i < end;i++) {
                localorder[i-i0] = i-i0;
            }
            vec_sort_by(localorder, sizeof(int64_t), crdiff+i0, end - i0 - 1);
            tmppivot = _find_breakpoint(crdiff + i0, end - i0 - 1);
            pivot = veci64_med(localorder + tmppivot, end - i0 - 1 - tmppivot) + i0 + 1;
            njumps = 1;
            for (i = tmppivot;i < end - i0 - 1;i++) {
                if (crdiff[i] > 2 * (currow[pivot + 1] - currow[pivot]))njumps++;
            }
            tmppivot = veci64_min(localorder + tmppivot, end - i0 - 1 - tmppivot) + i0 + 1;
            curminjump = currow[tmppivot + 1] - currow[tmppivot];
        }
        else {
            pivot = vec_maxid(crdiff + i0, end - i0 - 1) + i0 + 1;
        }
        
        if (steptype >= 1 && pivot>i0 + 2) {
            njumps=0;
            /*sum all internal differences for each row*/
            for(i=i0;i<pivot;i++){
                intradiff[i]=0;
                for(j=i0;j<pivot;j++){
                    intradiff[i] -= data[order[i]][order[j]] * data[order[i]][order[j]];
                }
            }
            /*find the "central" point and use it to extend the group*/
            tmppivot = pivot-1;
            vec_sort_by(order + i0, sizeof(int64_t), intradiff + i0, pivot - i0);
            for (i = pivot-1;i < end;i++) {
                currow[i] = data[order[pivot-1]][order[i]];
            }
            vec_sort_by(order + pivot - 1, sizeof(int64_t), currow + pivot - 1, end - pivot + 1);
            maxdiff=vec_max(currow+i0,pivot-i0);
            while(pivot < end && currow[pivot]<maxdiff)pivot++;
            steptype = 0;
            /* put the "central" point at the beginning of the cluster*/
            tmpi64 = order[tmppivot];
            order[tmppivot] = order[i0];
            order[i0] = tmpi64;
            steptype = 0;
            if (steptype == 2) {
                i0 = pivot;
                steptype--;
            }
        }
        else {
            steptype = 1;
            if (njumps >= 1 && !(pivot == end && i0 == begin)) {
                _distance_matrix_roughclust(data, i0, pivot, defaultcut, curminjump, &order);
            }
            if (end > pivot) {
                /*make sure the next group is close to the last */
                steptype = 2;
            }
            i0 = pivot;
            nclusters++;
        }
    }
    nclusters++;
    free(currow);
    free(crdiff);
    free(intradiff);
    free(localorder);

    *curorder = order;
    return  nclusters;
}
double _Ward_dist(double** data, size_t merged1, size_t merged2, size_t other, double w1, double w2, double wo) {
    double denom;
    double num;
    denom = w1 + w2 + wo;
    num = (w1+wo)*(data[merged1][other]) + (w2+wo)*(data[merged2][other]) - wo*(data[merged1][merged2]);
    return num / denom;
}
size_t _distance_matrix_wardclust(double** data, size_t begin, size_t end, double defaultcut, double minjump, int64_t** curorder, char** names, PF_t* output_tree) {
    size_t* parentclust;
    size_t* childrenclust;
    double** cdist;
    int64_t* order;
    size_t i, j, k;
    size_t id1, id2, idother;
    size_t size;
    size_t tmpsz;
    size_t to_merge;
    double* closestdist;
    double* F;
    double* D;
    double intravar, extravar;
    size_t* closestclust;
    size_t* curclustid;
    size_t* csize;
    size_t maxclustid;
    int has_unexplored_children;
    size_t nclust;
    double newdist;

    size = end - begin;

    order = *curorder;
    if (!order) {
        order = (int64_t*)malloc(sizeof(int64_t)*size);
        for (i = 0;i < size;i++) order[i] = (int64_t)(i);
    }
    parentclust = (size_t*)malloc(sizeof(size_t)*size * 2);
    childrenclust = (size_t*)malloc(sizeof(size_t)*size * 2);
    csize = (size_t*)calloc(size * 2, sizeof(size_t));
    F = (double*)malloc(sizeof(double)*size * 2);
    D = (double*)malloc(sizeof(double)*size * 2);
    curclustid = (size_t*)malloc(sizeof(size_t)*size);
    for (i = 0;i < size;i++) curclustid[i] = i;
    for (i = 0;i < size;i++) csize[i] = 1;

    cdist = (double**)malloc(sizeof(double*)*size);
    for (i = 0;i < size;i++) cdist[i] = (double*)malloc(sizeof(double)*size);
    extravar = 0;
    intravar = 0;
    for (i = 0;i < size;i++) {
        for (j = 0;j < size;j++) {
            cdist[i][j] = data[order[i + begin]][order[j + begin]] * data[order[i + begin]][order[j + begin]];
            extravar += cdist[i][j];
        }
    }
    extravar = extravar / 2;
    if (!(*curorder)) {
        order = (int64_t*)malloc(sizeof(int64_t)*size);
        for (i = 0;i < size;i++) order[i] = (int64_t)(i);
    }
    /* The rows and columns of the distance matrix are not shuffled throughout the procedure,
    the individual values are changed as clusters are merged together. 'order' is used to
    track which row/column corresponds to which of the actively modified clusters.
    In prallel, a dendrogram is created trhough parentid, and sizes of each cluster are
    stored in the same way. In order to navigate this data, 'clusterid' tracks which
    id is associated with the actively modified clusters.
    */
    closestdist = NULL;
    closestclust = NULL;
    closest_points(cdist, 0, size, order, &closestclust, &closestdist);
    nclust = 0;
    /* At the end of the following loop, we have produced a dendrogram stored in parentclust and childrenclust */
    for (i = 0;i < size - 1;i++) {
        to_merge = vec_minid(closestdist, size - i);
        id1 = to_merge;
        id2 = closestclust[to_merge];
        if (id2 < id1) {
            idother = id1;
            id1 = id2;
        }
        else {
            idother = id2;
        }
        D[curclustid[id1]] = sqrt(closestdist[to_merge]) * csize[curclustid[id1]] / (double)(csize[curclustid[id1]] + csize[curclustid[idother]]);
        D[curclustid[idother]] = sqrt(closestdist[to_merge]) * csize[curclustid[idother]] / (double)(csize[curclustid[id1]] + csize[curclustid[idother]]);
        intravar += (csize[curclustid[id1]] + csize[curclustid[idother]])*closestdist[to_merge];
        /* put id2 in last position and update distances */
        id2 = size - i - 1;
        tmpsz = order[idother]; order[idother] = order[id2]; order[id2] = tmpsz;
        tmpsz = curclustid[idother]; curclustid[idother] = curclustid[id2]; curclustid[id2] = tmpsz;
        closestclust[idother] = closestclust[id2];
        closestdist[idother] = closestdist[id2];

        closestdist[id1] = -1;
        extravar -= cdist[order[id1]][order[id2]];
        for (j = 0;j < id2;j++) {
            if (j != id1) {
                newdist = _Ward_dist(cdist, order[id1], order[id2], order[j], (double)csize[curclustid[id1]], (double)csize[curclustid[id2]], (double)csize[curclustid[j]]);
                extravar -= cdist[order[id1]][order[j]];
                extravar -= cdist[order[id2]][order[j]];
                extravar += newdist;
                if (closestdist[id1] == -1 || newdist < closestdist[id1]) {
                    closestdist[id1] = newdist;
                    closestclust[id1] = j;
                }
                cdist[order[id1]][order[j]] = newdist;
                cdist[order[j]][order[id1]] = newdist;
                if (closestdist[j] > newdist) {
                    closestdist[j] = newdist;
                    closestclust[j] = id1;
                }
                else if (closestclust[j] == id1 || closestclust[j] == idother || closestclust[j] == id2) {
                    closestclust[j] = id1;
                    closestdist[j] = newdist;
                    for (k = 0;k < id2;k++) {
                        if (k != j && closestdist[j] > cdist[order[j]][order[k]]) {
                            closestclust[j] = k;
                            closestdist[j] = cdist[order[j]][order[k]];
                        }
                    }
                }
            }
        }
        /* update modified clusters and their sizes */
        parentclust[curclustid[id1]] = size + i;
        parentclust[curclustid[id2]] = size + i;
        childrenclust[i * 2] = curclustid[id1];
        childrenclust[i * 2 + 1] = curclustid[id2];
        csize[size + i] = csize[curclustid[id1]] + csize[curclustid[id2]];
        curclustid[id1] = size + i;
        if (extravar > 0) {
            F[size + i] = intravar / extravar; /* the F statistic */
            /*F[size + i] = test_snedecor(intravar / extravar, size - i - 1, i); /* the associated p-value*/
            if (F[size + i] < defaultcut) nclust = size + i;
        }
    }
    maxclustid = size + i - 1;
    nclust = maxclustid - nclust + 1;

    /* do the sorting */
    k = maxclustid;
    while (k >= size) {
        id1 = childrenclust[(k - size) * 2];
        id2 = childrenclust[(k - size) * 2 + 1];
        if (id1 < id2) k = id1;
        else k = id2;
    }
    order[0] = (int64_t)k;
    for (i = 1; i < size;i++) {
        j = k;
        has_unexplored_children = 0;
        while (!has_unexplored_children && j != maxclustid) {
            k = j;
            j = parentclust[k];
            id1 = childrenclust[(j - size) * 2];
            id2 = childrenclust[(j - size) * 2 + 1];
            has_unexplored_children = (k < id1 || k < id2);
        }
        if (k < id1) k = id1;
        else if (k < id2) k = id2;
        while (k >= size) {
            id1 = childrenclust[(k - size) * 2];
            id2 = childrenclust[(k - size) * 2 + 1];
            if (id1 < id2) k = id1;
            else k = id2;
        }
        order[i] = (int64_t)k;
    }
    _print_dendrogram_as_newick_tree(parentclust, childrenclust, childrenclust + 1, D, F, names, size, output_tree, DENDRO_INTERLACED_CHILDREN);
    if (*curorder) {
        vec_reorder_byi64(*curorder, sizeof(int64_t), order, size);
        free(order);
    }
    else {
        *curorder = order;
    }

    for (i = 0;i < size;i++) free(cdist[i]);
    free(cdist);
    free(csize);
    free(parentclust);
    free(curclustid);
    free(closestclust);
    free(closestdist);
    free(childrenclust);
    free(F);
    free(D);
    return nclust;
}

size_t datatable_dmrclust(datatable_t* target, double defaultcut) {
    size_t res;
    size_t i;
    int64_t* neworder;
    if (!target || target->w != target->h) return 0;

    neworder = NULL;
    res = _distance_matrix_roughclust(target->values, 0, target->w, defaultcut, -1, &neworder);
    /*res = _distance_matrix_slcsort(target->values, target->w, &neworder);*/
    vec_reorder_byi64(target->colnames, sizeof(char*), neworder, target->w);
    vec_reorder_byi64(target->rownames, sizeof(char*), neworder, target->w);
    for (i = 0; i < target->w;i++) {
        DM64_assign(target->colids, target->colnames[i], (int)strlen(target->colnames[i]), (int64_t)i);
        DM64_assign(target->rowids, target->rownames[i], (int)strlen(target->rownames[i]), (int64_t)i);
        vec_reorder_byi64(target->values[i], sizeof(double), neworder, target->w);
    }
    vec_reorder_byi64(target->values, sizeof(double*), neworder, target->w);
    /*
    res = _dt_distance_matrix_roughsort(target, 0, target->w, target->w, defaultcut, -1);
    for (i = 0; i < target->w;i++) {
        DM64_assign(target->colids, target->colnames[i], (int)strlen(target->colnames[i]), (int64_t)i);
        DM64_assign(target->rowids, target->rownames[i], (int)strlen(target->rownames[i]), (int64_t)i);
    }
    */
    free(neworder);
    return res;
}
size_t datatable_dmwardclust(datatable_t* target, double defaultcut, PF_t* output_tree) {
    size_t res;
    size_t i;
    int64_t* neworder;
    if (!target || target->w != target->h) return 0;

    neworder = NULL;
    res = _distance_matrix_wardclust(target->values, 0, target->w, defaultcut, -1, &neworder, target->colnames, output_tree);
    vec_reorder_byi64(target->colnames, sizeof(char*), neworder, target->w);
    vec_reorder_byi64(target->rownames, sizeof(char*), neworder, target->w);
    for (i = 0; i < target->w;i++) {
        DM64_assign(target->colids, target->colnames[i], (int)strlen(target->colnames[i]), (int64_t)i);
        DM64_assign(target->rowids, target->rownames[i], (int)strlen(target->rownames[i]), (int64_t)i);
        vec_reorder_byi64(target->values[i], sizeof(double), neworder, target->w);
    }
    vec_reorder_byi64(target->values, sizeof(double*), neworder, target->w);
    free(neworder);
    return res;
}
size_t datatable_dmsingleclust(datatable_t* target, PF_t* output_tree) {
    size_t res;
    size_t i;
    int64_t* neworder;
    if (!target || target->w != target->h) return 0;

    neworder = NULL;
    res = _distance_matrix_slclust(target->values, target->w, &neworder, output_tree);
    vec_reorder_byi64(target->colnames, sizeof(char*), neworder, target->w);
    vec_reorder_byi64(target->rownames, sizeof(char*), neworder, target->w);
    for (i = 0; i < target->w;i++) {
        DM64_assign(target->colids, target->colnames[i], (int)strlen(target->colnames[i]), (int64_t)i);
        DM64_assign(target->rowids, target->rownames[i], (int)strlen(target->rownames[i]), (int64_t)i);
        vec_reorder_byi64(target->values[i], sizeof(double), neworder, target->w);
    }
    vec_reorder_byi64(target->values, sizeof(double*), neworder, target->w);
    free(neworder);
    return res;
}
size_t datatable_dmsortclust(datatable_t* target, double defaultcut) {
    size_t res;
    size_t i;
    int64_t* neworder;
    if (!target || target->w != target->h) return 0;

    neworder = NULL;
    res = _distance_matrix_progressivesort(target->values, 0, target->w, defaultcut, &neworder);
    vec_reorder_byi64(target->colnames, sizeof(char*), neworder, target->w);
    vec_reorder_byi64(target->rownames, sizeof(char*), neworder, target->w);
    for (i = 0; i < target->w;i++) {
        DM64_assign(target->colids, target->colnames[i], (int)strlen(target->colnames[i]), (int64_t)i);
        DM64_assign(target->rowids, target->rownames[i], (int)strlen(target->rownames[i]), (int64_t)i);
        vec_reorder_byi64(target->values[i], sizeof(double), neworder, target->w);
    }
    vec_reorder_byi64(target->values, sizeof(double*), neworder, target->w);
    free(neworder);
    return res;
}


void datatable_write(datatable_t* dt, PF_t* target, const char* sep, uint32_t flags) {
    size_t i, j;
    if (flags & DATATABLE_WRITECOLNAMES) {
        if (flags & DATATABLE_WRITEROWNAMES) {
            PFprintf(target, "%s", sep);
        }
        if (dt->w > 0) {
            for (j = 0;j < dt->w - 1;j++) {
                if(dt->colnames[j])
                    PFprintf(target, "%s", dt->colnames[j]);
                PFprintf(target, "%s", sep);
            }
            if (dt->colnames[j])
                PFprintf(target, "%s", dt->colnames[j]);
        }
        PFprintf(target, "\n");
    }
    for (i = 0;i < dt->h;i++) {
        if (flags & DATATABLE_WRITEROWNAMES) {
            if (dt->rownames[i])
                PFprintf(target, "%s", dt->rownames[i]);
            PFprintf(target, "%s", sep);
        }
        for (j = 0;j < dt->w - 1;j++) {
            for (j = 0;j < dt->w - 1;j++) {
                PFprintf(target, "%f%s", dt->values[i][j], sep);
            }
            if (dt->colnames[j])
                PFprintf(target, "%f", dt->values[i][j]);
        }
        PFprintf(target, "\n");
    }
}
