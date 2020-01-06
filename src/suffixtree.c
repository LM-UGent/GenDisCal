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
#include <stdlib.h>
#include <stdio.h>
#include "hash_table.h"
#include "suffixtree.h"
#include "vecops.h"

typedef struct bisize_htreenode {
    size_t v1;
    size_t v2;
    struct bisize_htreenode* parent;
    struct bisize_htreenode* backlink;
    struct bisize_htreenode** children;
    size_t nchild;
    size_t nchildalloc;
    size_t own_childid;
} bshtnode_t;

inline static bshtnode_t* bshtnode_alloc(void* datablock) {
    bshtnode_t* res;
    res = (bshtnode_t*) calloc(1, sizeof(bshtnode_t));
    return res;
}
inline static bshtnode_t* bshtnode_addnewchild(bshtnode_t* target, size_t v1, size_t v2) {
    if (target->nchildalloc < target->nchild+1) {
        if (target->nchildalloc == 0)target->nchildalloc = 2;
        else target->nchildalloc *= 2;
        target->children = (bshtnode_t**)realloc(target->children, sizeof(bshtnode_t*)*target->nchildalloc);
        memset(target->children+target->nchild,0,sizeof(bshtnode_t*)*(target->nchildalloc-target->nchild));
    }
    target->children[target->nchild] = bshtnode_alloc(NULL);
    target->children[target->nchild]->v1 = v1;
    target->children[target->nchild]->v2 = v2;
    target->children[target->nchild]->parent = target;
    target->children[target->nchild]->own_childid = target->nchild;
    target->nchild++;
    return target->children[target->nchild-1];
}
inline static void bshtnode_recountchildren(bshtnode_t* target) {
    size_t i, j;
    j = 0;
    for (i = 0;i < target->nchild;i++) {
        if (target->children[i]) {
            target->children[j] = target->children[i];
            target->children[j]->own_childid = j;
            j++;
        }
    }
    target->nchild = j;
}
static void bshtnode_addmissinglink(bshtnode_t* parent, bshtnode_t* child, size_t v1, size_t v2) {
    if (child->parent != parent) {
        if (child->parent) {
            child->parent->children[child->own_childid] = NULL;
            bshtnode_recountchildren(child->parent);
        }
        child->parent = bshtnode_addnewchild(parent, v1, v2);
    }
    else {
        parent->children[child->own_childid] = bshtnode_alloc(NULL);
        parent->children[child->own_childid]->v1 = v1;
        parent->children[child->own_childid]->v2 = v2;
        parent->children[child->own_childid]->parent = parent;
        parent->children[child->own_childid]->own_childid = child->own_childid;
        child->parent = parent->children[child->own_childid];
    }
    child->parent->children = (bshtnode_t**)malloc(sizeof(bshtnode_t*));
    child->parent->children[0] = child;
    child->parent->nchild = 1;
    child->parent->nchildalloc = 1;
    child->own_childid = 0;
}
inline static void bshtnode_free(bshtnode_t* target, int recalculate_parent) {
    size_t i;
    bshtnode_t* caller;
    size_t parentchildren;
    caller = target;
    if (target->parent){
        if (recalculate_parent == 1) {
            target = target->parent;
            parentchildren = target->nchild;
            if (caller->nchild > 0) {
                while (target->nchildalloc < target->nchild + caller->nchild) {
                    if (target->nchildalloc == 0)target->nchildalloc = 2;
                    else target->nchildalloc *= 2;
                    target->children = (bshtnode_t**)realloc(target->children, sizeof(bshtnode_t*)*target->nchildalloc);
                    memset(target->children + target->nchild, 0, sizeof(bshtnode_t*)*(target->nchildalloc - target->nchild));
                }
                target->children[caller->own_childid] = NULL;
                bshtnode_recountchildren(target);
                parentchildren = target->nchild;
                memcpy(target->children + target->nchild, caller->children, sizeof(bshtnode_t*)*caller->nchild);
                target->nchild += caller->nchild;
            }
        }
        else {
            target->parent->children[target->own_childid] = NULL;
        }
    }
    else {
        target = NULL;
        parentchildren = 0;
    }
    /* note : if at this point target==NULL and caller->nchild > 0, the following may result in 'lost' children
              (i.e. memory leak). However, if the nodes are reference elsewhere this may not be case, so the program
              will not be terminated here.
    */
    if (caller->nchild > 0) {
        for (i = 0;i < caller->nchild;i++) {
            if (caller->children[i]) {
                caller->children[i]->parent = target;
                caller->children[i]->own_childid = parentchildren + i;
            }
        }
        free(caller->children);
    }
    free(caller);
}
static void bshtnode_cutbranch(bshtnode_t* target, int affects_parent) {
    size_t depth;
    size_t i;
    depth = 0;
    while (target->nchild == 1) {
        depth++;
        target = target->children[0];
    }
    for (i = 0;i < target->nchild;i++) {
        if (target->children[i])
            bshtnode_cutbranch(target->children[i], 0);
    }
    while (depth > 0) {
        target = target->parent;
        bshtnode_free(target->children[0], 0);
        depth--;
    }
    bshtnode_free(target, affects_parent);
}

struct suffixtree_t {
    bshtnode_t* root;
    char* data;
    size_t datalen;
    int freedata, sorted;
    size_t longest_repeat;
};


void fprintmarks(FILE* fp, char* marks, size_t datalen) {
    size_t i;
    for (i = 0;i < datalen;i++) {
        if (marks[i])fputc('O', fp);
        else fputc('-', fp);
    }
    fputc('\n', fp);
}


static int suftree_printnode(bshtnode_t* target, suftree_t* tree,int format) {
    size_t i;
    size_t indent;
    indent = 0;
    if (format % 2 == 0) {
        if (target->v2 - target->v1 < 10) {
            for (i = target->v1;i < target->v2;i++) {
                putchar(tree->data[i]);
            }
            indent = target->v2 - target->v1;
        }
        else {
            for (i = 0; i < 4;i++) { putchar(tree->data[i + target->v1]); }
            printf("..");
            for (i = 0; i < 4;i++) { putchar(tree->data[i + target->v2 - 4]); }
            indent = 10;
        }
        if (target->v2 == tree->datalen || target->nchild == 0) {
            indent++;
            putchar('$');
            if (format == 0)
                putchar('\n');
        }
    }
    if (format == 1) {
        for (i = target->v1;i < target->v2;i++) {
            putchar(tree->data[i]);
            indent = target->v2 - target->v1;
        }
        indent = target->v2 - target->v1;
        if (target->v2 == tree->datalen || target->nchild == 0) {
            indent++;
            putchar('$');
        }
    }
    return (int) indent;
}
static void suftree_printprefix(bshtnode_t* target, suftree_t* tree) {
    if (target->parent) {
        suftree_printprefix(target->parent,tree);
    }
    suftree_printnode(target,tree, 0);
}
static int suftree_printnode_and_backlink(bshtnode_t* target, suftree_t* tree) {
    size_t i;
    size_t indent;
    indent = 0;
    if (target->parent)
        suftree_printprefix(target->parent,tree);
    if (target->v2 - target->v1 < 10) {
        for (i = target->v1;i < target->v2;i++) {
            putchar(tree->data[i]);
        }
        indent = target->v2 - target->v1;
    }
    else {
        for (i = 0; i < 4;i++) { putchar(tree->data[i + target->v1]); }
        printf("..");
        for (i = 0; i < 4;i++) { putchar(tree->data[i + target->v2 - 4]); }
        indent = 10;
    }
    if (target->v2 == tree->datalen || target->nchild == 0) {
        indent++;
        putchar('$');
    }
    printf("@" _LLD_, target->v2);
    printf(" -> ");
    if (target->backlink) {
        target = target->backlink;
        if (target->parent)
            suftree_printprefix(target->parent, tree);
        if (target->v2 - target->v1 < 10) {
            for (i = target->v1;i < target->v2;i++) {
                putchar(tree->data[i]);
            }
            indent = target->v2 - target->v1;
        }
        else {
            for (i = 0; i < 4;i++) { putchar(tree->data[i + target->v1]); }
            printf("..");
            for (i = 0; i < 4;i++) { putchar(tree->data[i + target->v2 - 4]); }
            indent = 10;
        }
        if (target->v2 == tree->datalen || target->nchild == 0) {
            indent++;
            putchar('$');
        }
        printf("@" _LLD_, target->v2);
    }
    putchar('\n');
    return (int) indent;
}
void suftree_printbacklinks(suftree_t* tree) {
    int morenodes;
    size_t curchild;
    bshtnode_t* curnode;
    morenodes = 0;
    curnode = tree->root;
    if (curnode->nchild > 0) {
        morenodes = 1;
        curnode = tree->root->children[0];
    }
    curchild = 0;
    while (morenodes) {
        if(curchild == 0 )
            suftree_printnode_and_backlink(curnode, tree);
        if (curchild < curnode->nchild) {
            curnode = curnode->children[curchild];
            curchild = 0;
        }
        else {
            if (curnode->parent) {
                curchild = curnode->own_childid + 1;
                curnode = curnode->parent;
            }
            else
                morenodes = 0;
        }
    }
}
char* suftree_getstrptr(suftree_t* tree, size_t* output_len) {
    if (!tree) {
        *output_len = 0;
        return NULL;
    }
    if (output_len) *output_len = tree->datalen;
    return tree->data;
}
size_t suftree_count_node_descendants(suftree_t* tree, bshtnode_t* node) {
    size_t result;
    size_t i;
    result = node->nchild;
    for (i = 0;i < node->nchild;i++)
        result += suftree_count_node_descendants(tree, node->children[i]);
    return result;
}
void suftree_node_addresses_to_ids(ht64_t* data, bshtnode_t* curnode, int64_t* curid) {
    int nf;
    size_t i;
    ht64_set(data, (char*)(&curnode), sizeof(bshtnode_t*), *curid, &nf);
    *curid = (*curid) + 1;
    for (i = 0;i < curnode->nchild;i++)
        suftree_node_addresses_to_ids(data, curnode->children[i], curid);
}
void suftree_write_node_to_file(bshtnode_t* node, PF_t* file, ht64_t* addr2id) {
    int nf;
    int64_t id;
    PFputint64(file, node->v1);
    PFputint64(file, node->v2);
    PFputint64(file, node->nchild);
    id = ht64_get(addr2id, (char*)(&(node->backlink)), sizeof(bshtnode_t*), &nf);
    if (nf) id = 0;
    else id += 1;
    PFputint64(file, id);
}
void suftree_save(suftree_t* tree, PF_t* target)
{
    int stop;
    size_t numnodes;
    bshtnode_t* curnode;
    ht64_t* addr2id;
    int64_t curid;
    PFputint64(target,tree->datalen);
    PFputint16(target, (int16_t)tree->sorted);
    PFwrite(tree->data, 1, tree->datalen, target);
    stop = 0;
    numnodes = 0;
    curnode = tree->root;
    numnodes = suftree_count_node_descendants(tree, curnode) + 1;
    PFputint64(target, (int64_t) numnodes);
    addr2id = ht64_alloc_size(numnodes);
    curid = 0;
    suftree_node_addresses_to_ids(addr2id, curnode, &curid);
    
    ht64_free(addr2id);
}
void suftree_parse_loaded_nodes(suftree_t* target, bshtnode_t** loadednodes, size_t* backlinks, size_t* curnodeid, bshtnode_t* parentnode, size_t childid) {
    size_t i;
    size_t nodeid;
    nodeid = *curnodeid;
    if (parentnode == NULL) {
        target->root = loadednodes[nodeid];
    }
    else {
        parentnode->children[childid] = loadednodes[nodeid];
    }
    loadednodes[nodeid]->parent = parentnode;
    loadednodes[nodeid]->children = calloc(loadednodes[nodeid]->nchild, sizeof(bshtnode_t*));
    loadednodes[nodeid]->nchildalloc = loadednodes[nodeid]->nchild;
    if (backlinks[nodeid] > 0)
        loadednodes[nodeid]->backlink = loadednodes[backlinks[nodeid] - 1];
    *curnodeid = (*curnodeid) + 1;
    for (i = 0;i < loadednodes[nodeid]->nchild;i++) {
        suftree_parse_loaded_nodes(target, loadednodes, backlinks, curnodeid, loadednodes[nodeid], i);
    }
}
suftree_t* suftree_load(PF_t* target)
{
    suftree_t* result;
    bshtnode_t** id2addr;
    size_t* backlinks;
    size_t i,numnodes;
    result = (suftree_t*) calloc(1, sizeof(suftree_t));
    result->datalen = PFgetint64(target);
    result->sorted = PFgetint16(target);
    result->data = malloc(result->datalen);
    PFread(result->data, 1, result->datalen, target);
    numnodes = (size_t)PFgetint64(target);
    id2addr = (bshtnode_t**)calloc(numnodes, sizeof(bshtnode_t*));
    backlinks = (size_t*)calloc(numnodes, sizeof(size_t));
    for (i = 0;i < numnodes;i++) {
        id2addr[i] = bshtnode_alloc(NULL);
        id2addr[i]->v1 = (size_t)PFgetint64(target);
        id2addr[i]->v2 = (size_t)PFgetint64(target);
        id2addr[i]->nchild = (size_t)PFgetint64(target);
        backlinks[i] = (size_t)PFgetint64(target);
    }
    i = 0;
    suftree_parse_loaded_nodes(result, id2addr, backlinks, &i, NULL, 0);
    result->freedata = 1;
    return result;
}
static void suftree_printbranch(bshtnode_t* target, int indent, suftree_t* tree, int format) {
    size_t i;
    int indent0;
    int j;
    bshtnode_t* startnode;
    if (format == 0) {
        for (j = 0;j < indent;j++) {
            putchar(' ');
        }
        if (indent < 0)indent = -indent;
        indent += suftree_printnode(target, tree, 0);
        while (target->nchild == 1) {
            target = target->children[0];
            suftree_printnode(target, tree, 0);
            indent += (int)(target->v2 - target->v1);
        }
        for (i = 0;i < target->nchild;i++) {
            suftree_printbranch(target->children[i], indent*(i ? 1 : -1), tree, 0);
        }
    }
    if (format == 1) {
        startnode = target;
        if (indent < 0)indent = -indent;
        indent0 = indent;
        if (target->nchild > 0) {
            indent += (int)(target->v2 - target->v1);
            while (target->nchild == 1) {
                target = target->children[0];
                indent += (int)(target->v2 - target->v1);
            }
            putchar('(');
            for (i = 0;i < target->nchild;i++) {
                if (i > 0)putchar(',');
                suftree_printbranch(target->children[i], indent*(i ? 1 : -1), tree, 1);
            }
            putchar(')');
            target = startnode;
            indent = indent0;
            while (target->nchild == 1) {
                target = target->children[0];
                suftree_printnode(target, tree, 1);
                printf("@%d", indent);
                indent += (int)(target->v2 - target->v1);
            }
        }
        suftree_printnode(target, tree, 1);
        printf("@%d", indent);
    }
}
void suftree_printf(suftree_t* target) {
    suftree_printbranch(target->root, 0, target, 0);
}
void suftree_printnewick(suftree_t* target) {
    suftree_printbranch(target->root, 0, target, 1);
}

static char suftreenode_firstletter(suftree_t* tree, bshtnode_t* node, int* nullflag) {
    if (node->v1 >= tree->datalen || node->v2 == 0) {
        *nullflag = 1;
        return 0;
    }
    return tree->data[node->v1];
}
static void suftreenode_sort(suftree_t* tree, bshtnode_t* target) {
    size_t i,j, maxi;
    size_t targi;
    bshtnode_t* tmp;
    maxi = target->nchild;
    /* have 0 or 1 child means that the node is inherently sorted, so do nothing */
    if (maxi < 2) return;
    /*
    insertion sort should work fine here considering we expect a low number of children per node
    At worst, we could expect n=257, but such nodes will not be common
    */
    while (maxi>0 && target->children[maxi - 1]->v1 >= tree->datalen || target->children[maxi - 1]->v2 == 0) {
        maxi--;
    }
    i = 0;
    while (i<maxi) {
        if (target->children[i]->v1 >= tree->datalen || target->children[i]->v2 == 0) {
            targi = maxi;
            tmp = target->children[maxi - 1];
            target->children[maxi - 1] = target->children[i];
            target->children[i] = tmp;
            while (maxi>0 && target->children[maxi - 1]->v1 >= tree->datalen || target->children[maxi - 1]->v2 == 0) {
                maxi--;
            }
        }
        else {
            targi = i;
            while (targi > 0 && tree->data[target->children[i]->v1] < tree->data[target->children[targi-1]->v1]) {
                targi--;
            }
            tmp = target->children[i];
            for (j = i;j > targi; j--) {
                target->children[j] = target->children[j - 1];
            }
            target->children[targi] = tmp;
            i++;
        }
    }
    for (i = 0;i < target->nchild;i++) {
        target->children[i]->own_childid = i;
    }
}
static void suftreenode_recursivesort(suftree_t* tree, bshtnode_t* target) {
    size_t i;
    suftreenode_sort(tree, target);
    for (i = 0;i < target->nchild;i++) {
        suftreenode_recursivesort(tree, target->children[i]);
    }
}
void suftree_sort(suftree_t* tree) {
    suftreenode_recursivesort(tree, tree->root);
    tree->sorted = 1;
}
size_t suftree_nodeprefixlen(bshtnode_t* node, suftree_t* tree) {
    size_t result;
    bshtnode_t* curnode = node->parent;
    result = 0;
    while (curnode) {
        result += curnode->v2 - curnode->v1;
        curnode = curnode->parent;
    }
    return result;
}
static inline void suftree_naive_build(suftree_t* tree) {
    /*
        The naive algorithm works pretty well until you begin to have big repeats in you input, since the cost of each loop is O(m)
        where m is the length of the longest previously encountered repeated prefix. Under random conditions, m ~ log(a,n), where a is the alphabet size
    */
    size_t i, j, k, m;
    size_t previnternodedepth;
    size_t curnodedepth;
    size_t tmpnodedepth;
    bshtnode_t* curnode;
    bshtnode_t* prevnode;
    bshtnode_t* tmpnode;
    bshtnode_t* previnternode;
    int newchild;
    int splitnode;
    int omega = 0;
    size_t nnodes;
    nnodes = 0;
    tree->root = bshtnode_alloc(NULL);
    prevnode = NULL;
    tmpnode = NULL;
    previnternode = NULL;
    for (i = 0;i < tree->datalen;i++) {
        m = i;
        newchild = 0;
        curnode = tree->root;
        j = k = splitnode = 0;
        curnodedepth = 0;
        while (!newchild)
        {
            splitnode = 0;
            newchild = 1;
            for (j = 0;j < curnode->nchild; j++) {
                k = curnode->children[j]->v1;
                if (tree->data[k] == tree->data[m]) {
                    newchild = 0;
                    break;
                }
            }
            curnodedepth += curnode->v2 - curnode->v1;
            if (!newchild) {
                splitnode = 1;
                curnode = curnode->children[j];
                while (k < curnode->v2) {
                    if (m>=tree->datalen || tree->data[k] != tree->data[m]) {
                        curnodedepth += k - curnode->v1;
                        newchild = 1;
                        break;
                    }
                    k++;
                    m++;
                }
            }
        }
        /* each pass should result in exactly 1 new child (with potentially one intermediary node */
        if (splitnode) {
            /* 
            the creation of the intermediary node must be done first,
            since the next step might create a backlink to it
            */
            bshtnode_addmissinglink(curnode->parent, curnode, curnode->v1, k);
            tmpnode = curnode->parent;
            curnode->v1 = k;
            curnode = tmpnode;
            nnodes++;
        }
        if (prevnode) {
            /* this only creates back-links for the tree leaves */
            prevnode->backlink = bshtnode_addnewchild(curnode, m, tree->datalen);
            prevnode = prevnode->backlink;
        }
        else
            prevnode = bshtnode_addnewchild(curnode, m, tree->datalen);
        if (previnternode && (previnternode->parent!=tree->root || previnternode->v2- previnternode->v1 > 1)) {
            tmpnode = curnode;
            tmpnodedepth = curnodedepth;
            while (tmpnode->parent && tmpnodedepth >= previnternodedepth) {
                tmpnodedepth -= tmpnode->v2 - tmpnode->v1;
                tmpnode = tmpnode->parent;
            }
            if (tmpnode != tree->root)
                previnternode->backlink = tmpnode;
            previnternode = NULL;
        }
        if (splitnode) {
            previnternode = curnode;
            previnternodedepth = curnodedepth;
        }
        nnodes++;
    }
}
suftree_t* suftree_from(char* data, int option, size_t datalen) {
    suftree_t* result;

    if ((option & 0b001111) == 0)
        option |= SUFTREE_NAIVE;
    if ((option & 0b110000) == 0)
        option |= SUFTREECOPYTEXT;
    
    result = (suftree_t*) calloc(1,sizeof(suftree_t));
    result->root = NULL;
    result->datalen = datalen;
    if (option & SUFTREECOPYTEXT) {
        result->data = (char*)malloc(result->datalen);
        memcpy(result->data, data, datalen);
        result->freedata = 1;
    }
    else {
        result->data = data;
        result->freedata = 0;
    }
    if ((option & 0b1111) == SUFTREE_NAIVE) {
        suftree_naive_build(result);
    }
    return result;
}
void suftree_free(suftree_t* target) {
    if (target->freedata && target->data) {
        free(target->data);
    }
    target->data = NULL;
    bshtnode_cutbranch(target->root, 0);
}

suftree_t* suftree_intersection(suftree_t* A, suftree_t* B) {
    suftree_t* result;
    bshtnode_t* curnodeA;
    bshtnode_t* curnodeB;
    bshtnode_t* curnodeR;
    bshtnode_t* relevantparentA;
    bshtnode_t* relevantparentB;
    size_t iA, iB, jA, jB, j, advanced;
    char cA, cB;
    int finished, mismatch,split;
    size_t depth;
    int nf;
    if (!A->sorted)
        suftree_sort(A);
    if (!B->sorted)
        suftree_sort(B);
    result = calloc(1, sizeof(suftree_t));
    result->data = A->data;
    result->datalen = A->datalen;
    result->root = bshtnode_alloc(NULL);
    
    curnodeR = result->root;
    curnodeA = A->root;
    curnodeB = B->root;
    finished = 0;
    iA = 0;
    iB = 0;
    depth = 0;
    while (!finished) {
        nf = 1;
        if (iA < curnodeA->nchild && iB < curnodeB->nchild) {
            cA = suftreenode_firstletter(A, curnodeA->children[iA], &nf);
            cB = suftreenode_firstletter(B, curnodeB->children[iB], &nf);
            nf = 0;
        }
        /* find the next matching child */
        while (iA < curnodeA->nchild && iB < curnodeB->nchild && cA != cB && nf==0) {
            while (iA < curnodeA->nchild && cA < cB && nf==0) {
                iA++;
                if (iA >= curnodeA->nchild)
                    nf = 1;
                else
                    cA = suftreenode_firstletter(A, curnodeA->children[iA], &nf);
            }
            while (iB < curnodeB->nchild && cB < cA && nf == 0) {
                iB++;
                if (iB >= curnodeB->nchild)
                    nf = 1;
                else
                    cB = suftreenode_firstletter(B, curnodeB->children[iB], &nf);
            }
        }
        /*
        at this point, either the iA'th child of A and iB'th child of B start with the first letter,
        or the analysis of this node is over
        */
        if (cA == cB && !nf) {
            relevantparentA = curnodeA;
            relevantparentB = curnodeB;
            curnodeA = curnodeA->children[iA];
            curnodeB = curnodeB->children[iB];
            advanced = 0;
            jA = curnodeA->v1;
            jB = curnodeB->v1;
            curnodeR = bshtnode_addnewchild(curnodeR, jA, jA);
            mismatch=split = 0;
            while ( !(mismatch || split)) {
                if (jA == curnodeA->v2) {
                    if (curnodeA->nchild == 1) {
                        if (curnodeA->children[0]->v2 != 0) {
                            curnodeA = curnodeA->children[0];
                            curnodeR->v1 = curnodeA->v1 - advanced;
                        }
                        else
                            mismatch = 1;
                    }
                    else if(jB < curnodeB->v2) {
                        cB = B->data[jB];
                        for (j = 0;j < curnodeA->nchild;j++) {
                            if (suftreenode_firstletter(A, curnodeA->children[j], &nf) == cB)
                                break;
                        }
                        if (j < curnodeA->nchild) {
                            curnodeA = curnodeA->children[j];
                            jA = curnodeA->v1;
                            curnodeR->v1 = curnodeA->v1 - advanced;
                        }
                        else {
                            mismatch = 1;
                        }
                    }
                    else {
                        split = 1;
                    }
                }
                else if (jB == curnodeB->v2) {
                    cA = A->data[jA];
                    if (curnodeB->nchild == 1) {
                        if (curnodeB->children[0]->v2 != 0)
                            curnodeB = curnodeB->children[0];
                        else
                            mismatch = 1;
                    }
                    else {
                        for (j = 0;j < curnodeB->nchild;j++) {
                            if (suftreenode_firstletter(B, curnodeB->children[j], &nf) == cA)
                                break;
                        }
                        if (j < curnodeB->nchild) {
                            curnodeB = curnodeB->children[j];
                            jB = curnodeB->v1;
                        }
                        else {
                            mismatch = 1;
                        }
                    }
                }
                if (!split && !mismatch) {
                    if (jB < B->datalen)
                        cB = B->data[jB];
                    if (jA < A->datalen)
                        cA = A->data[jA];
                }
                if (cA != cB)mismatch = 1;
                if (!(split || mismatch)) {
                    if (A->data[jA] != B->data[jB])mismatch = 1;
                    jA++;
                    jB++;
                    advanced++;
                    curnodeR->v2 = jA;
                    depth++;
                }
            }
            if (mismatch) {
                while (curnodeA->parent != relevantparentA)
                    curnodeA = curnodeA->parent;
                iA = curnodeA->own_childid + 1;
                curnodeA = curnodeA->parent;
                while (curnodeB->parent != relevantparentB)
                    curnodeB = curnodeB->parent;
                curnodeR = curnodeR->parent;
                iB = curnodeB->own_childid + 1;
                curnodeB = curnodeB->parent;
            }
            if (split) {
                iA = 0;
                iB = 0;
            }
        }
        else {
            if (result->longest_repeat < depth)
                result->longest_repeat = depth;
            iA = curnodeA->own_childid;
            jA = curnodeA->v2 - curnodeA->v1;
            curnodeA = curnodeA->parent;
            iB = curnodeB->own_childid;
            jB = curnodeB->v2 - curnodeB->v1;
            curnodeB = curnodeB->parent;
            if (jA == jB) {
                curnodeR = curnodeR->parent;
                iA++;
                iB++;
            }
            if (!curnodeA || !curnodeB)finished = 1;
            while ((jA!=jB || jA==0 || iA>=curnodeA->nchild || iB>=curnodeB->nchild) && finished == 0) {
                if (curnodeA->parent == NULL)finished = 1;
                else {
                    iA = curnodeA->own_childid;
                    jA += curnodeA->parent->v2 - curnodeA->parent->v1;
                    curnodeA = curnodeA->parent;
                    iB = curnodeB->own_childid;
                    jB += curnodeB->parent->v2 - curnodeB->parent->v1;
                    curnodeB = curnodeB->parent;
                    if (jA == jB) {
                        curnodeR = curnodeR->parent;
                        iA++;
                        iB++;
                    }
                }
            }
            depth -= jA;
        }
    }
    return result;
}
size_t suftree_search(suftree_t* target, char* data, size_t datalen) {
    size_t i,j,k;
    bshtnode_t* curnode;
    curnode = target->root;
    i = 0;
    for (j = 0;j < curnode->nchild;j++) {
        if (data[i] == target->data[curnode->children[j]->v1])
            break;
    }
    if (j == curnode->nchild)
        curnode = NULL;
    else {
        curnode = curnode->children[j];
        k = curnode->v1;
    }
    while (curnode && i < datalen && data[i] == target->data[k]) {
        while (k<curnode->v2 && i<datalen && data[i] == target->data[k]) {
            i++;
            k++;
        }
        if (k == curnode->v2 && i<datalen) {
            for (j = 0;j < curnode->nchild;j++) {
                if (data[i] == target->data[curnode->children[j]->v1])
                    break;
            }
            if (j == curnode->nchild)
                curnode = NULL;
            else {
                curnode = curnode->children[j];
                k = curnode->v1;
            }
        }
    }
    if(!curnode || i<datalen) return target->datalen;
    return k - datalen;
}
char* suftree_markpresent(char* data, size_t datalen, suftree_t* tree, size_t minlen, size_t* correspondance) {
    size_t i, j, k, matches, lastmatches, i0, lasti0;
    char* result;
    int mismatch;
    bshtnode_t* curnode;
    int new_good_fragment;
    
    result = calloc(datalen,1);
    i = 0;
    if (correspondance) {
        for (j = 0;j < datalen;j++) {
            correspondance[j] = tree->datalen;
        }
    }
    lastmatches = 0;
    lasti0 = 0;
    while (i < datalen) {
        /* find the first mismatch */
        curnode = tree->root;
        mismatch = 0;
        matches = 0;
        /*note: k = first mismatched character */
        k = 0;
        i0 = i;
        while (!mismatch) {
            mismatch = 1;
            if (k >= curnode->v2) {
                for (j = 0;j < curnode->nchild;j++) {
                    if (tree->data[curnode->children[j]->v1] == data[i]) {
                        mismatch = 0;
                        k = curnode->children[j]->v1;
                        curnode = curnode->children[j];
                        break;
                    }
                }
            }
            else {
                if (tree->data[k] == data[i]) {
                    mismatch = 0;
                    k++;
                    i++;
                    matches++;
                }
            }
            if (i == datalen) mismatch = 1;
        }
        /* mark result if it is within the minimum length */
        new_good_fragment = 0;
        if (matches > minlen) {
            if (i0 < lasti0 + lastmatches) {
                if (lastmatches > matches) matches -= (lasti0 + lastmatches - i0);
                if (matches > minlen) new_good_fragment = 1;
            }
            else {
                new_good_fragment = 1;
            }
            if (new_good_fragment) {
                memset(result + i - matches, 1, matches);
                if (i0 - lasti0 < minlen) {
                    for (j = 0;j < i0 - lasti0;j++) {
                        if (result[i0 - j - 1]) {
                            result[i0 - j - 1] = 0;
                            correspondance[i0 - j - 1] = tree->datalen;
                        }
                    }
                }
                if (correspondance) {
                    for (j = 0;j < matches;j++) {
                        correspondance[i - matches + j] = k - matches + j;
                    }
                }
                lasti0 = i0;
                lastmatches = matches;
            }
        }
        else if (i0 > lasti0 + lastmatches) {
            lastmatches = 0;
        }
        /* follow backlinks until reaching a point where the current position is no longer part of the branch */
        if (k > 0) {
            k--;
            if (curnode->v2 <= curnode->v1 && curnode->v2 > 0)
                curnode = curnode->parent;
            while (k >= curnode->v1 && curnode->backlink && curnode->backlink!=curnode) {
                curnode = curnode->backlink;
                i0++;
            }
            if (matches <= 1) i0++;
            i = i0;
        }
        else {
            i = (++i0);
        }
    }
    return result;
}
void dilate_1D(char* marks, size_t datalen, size_t amount, char color) {
    size_t i, j;
    i = 1;
    while (i < datalen) {
        if (marks[i] && !marks[i - 1]) {
            j = 0;
            while (j<amount && i - j>0) {
                j++;
                if (!marks[i - j])
                    marks[i - j] = color;
            }
        }
        if (!marks[i] && marks[i - 1]) {
            j = 0;
            while (j<amount && i + j < datalen) {
                if (!marks[i + j])
                    marks[i + j] = color;
                j++;
            }
            i += amount;
        }
        i++;
    }
}
void erode_1D(char* marks, size_t datalen, size_t amount, char color) {
    size_t i, j;
    i = 1;
    while (i < datalen-1) {
        if (!marks[i] && marks[i-1]) {
            j = 0;
            while (j<amount && i - j>0 && marks[i - j - 1]) {
                j++;
                marks[i - j] = color;
            }
        }
        if (marks[i] && !marks[i-1]) {
            j = 0;
            while (j<amount && i + j < datalen && marks[i + j]) {
                marks[i + j] = color;
                j++;
            }
            i += j;
        }
        i++;
    }
}
void bridgegaps_1D(char* marks, size_t datalen, size_t amount, char color) {
    size_t i, count;
    i = 1;
    count = 0;
    while (i < datalen - 1) {
        if (!marks[i] && marks[i - 1]) {
            count = 0;
        }
        if (marks[i] && !marks[i - 1]) {
            if (count <= amount) {
                memset(marks + i - count, color, count);
            }
        }
        if (!marks[i]) count++;
        i++;
    }
}

void suftree_printmarkeddata(char* marks, char* data, size_t datalen, suftree_t* tree, int outputformat) {
    size_t i, j, k;
    i = 1;
    if (outputformat == DATAMARK_BOTTOMSTARS) {
        fwrite(data, 1, datalen, stdout);
        for (i = 0;i < datalen;i++) {
            if (marks[i])marks[i] = '*';
        }
        printf("\n");
        fwrite(marks, 1, datalen, stdout);
    }
    else if (outputformat == DATAMARK_INTERVALS) {
        j = 0;
        k = 0;
        for (i = 0;i < datalen;i++) {
            if (!marks[i]) {
                if (i - j > 1) {
                    fwrite(data + j + 1, 1, i - j, stdout);
                }
                j = i;
            }
            else {
                if (j == i - 1) {
                    fprintf(stdout, "-[" _LLD_ "-" _LLD_ "]-", k, i);
                }
                k = i + 1;
            }
        }
        if (i - j > 1) {
            fwrite(data + j + 1, 1, i - j, stdout);
        }
        if (k < datalen) {
            fprintf(stdout, "-[" _LLD_ "-" _LLD_ "]-", k, i);
        }
    }
    printf("\n");
}

static inline int acceptablegap(size_t* data, size_t left, size_t right, size_t maxgaperror, size_t alignlen) {
    if (right - left > alignlen)
        return 0;
    if (data[left] < data[right])
        return ((data[right] - data[left]) < maxgaperror + (right - left));
    else
        return ((data[left] - data[right]) < maxgaperror + (right - left));
}

void bridge_substitutions(char* marks, size_t* correspondance, char* data, size_t datalen, char* source, size_t sourcelen, size_t minlen, double wiggle_room, size_t maxgap) {
    size_t pos, i;
    size_t readupto;
    size_t srcmove, qrymove, delta, avgmove;
    size_t lastblockend, lastblockstart, curblockstart;
    size_t matching_lbe, matching_lbs, matching_cbs;
    size_t lastgoodblockend;
    double slope;
    int is_good_block;
    double relative_advance;
    double expected_delta;
    double extension_tol;
    
    matching_lbe = 0;
    matching_lbs = 0;
    matching_cbs = 0;
    lastblockend = 0;
    lastblockstart = 0;
    curblockstart = 0;

    if (datalen < 1) {
        /*nothing happens*/
        return;
    }
    pos = 0;
    readupto = datalen - 1;
    while (pos < readupto && !marks[pos]) {
        pos++;
    }
    if (pos < readupto) {
        lastblockstart = pos;
        matching_lbs = correspondance[curblockstart];
    }
    /* generate good blocks */
    while (pos < readupto) {
        if (marks[pos] && !marks[pos + 1]) {
            if (pos + 1 - lastblockstart > minlen)
                lastgoodblockend = lastblockend;
            lastblockend = pos+1;
            matching_lbe = correspondance[pos]+1;
        }
        else if (!marks[pos] && marks[pos + 1]) {
            curblockstart = pos+1;
            matching_cbs = correspondance[curblockstart];
            qrymove = curblockstart-lastblockend;
            is_good_block = 0;
            if (matching_cbs > matching_lbe) {
                srcmove = matching_cbs - matching_lbe;
                avgmove = (srcmove + qrymove + 1) / 2;
                if (srcmove >= qrymove) delta = srcmove - qrymove;
                else delta = qrymove - srcmove;
                if ((double)delta < wiggle_room*(double)avgmove && avgmove<maxgap) {
                    for (i = lastblockend;i < curblockstart;i++) {
                        relative_advance = ((double)(i - lastblockend)) / ((double)qrymove);
                        marks[i] = 1;
                        correspondance[i] = matching_lbe + (size_t)(relative_advance*(double)srcmove);
                    }
                    is_good_block = 1;
                    curblockstart = lastblockstart;
                }
            }
            if (lastblockend - lastblockstart >= minlen) {
                is_good_block = 1;
            }
            if (!is_good_block) {
                if (lastblockend - lastblockstart < minlen) {
                    i = lastblockstart;
                    while (i < lastblockend) {
                        marks[i] = 0;
                        /*correspondance[i] = sourcelen;*/
                        i++;
                    }
                }
            }
            lastblockstart = curblockstart;
            matching_lbs = matching_cbs;
        }
        else if (marks[pos] && marks[pos + 1]) {
            if (pos + 1 < readupto) {
                if (correspondance[pos] >= correspondance[pos + 1])
                    delta = correspondance[pos] - correspondance[pos + 1] + 1;
                else
                    delta = correspondance[pos + 1] - correspondance[pos] - 1;
                if (delta > minlen) {
                    marks[pos + 1] = 0;
                    pos--;
                }
            }
        }
        pos++;
    }
    /* compare correspondance to expected value based on good blocks */
    pos = 0;
    while (pos < readupto) {
        if (!marks[pos] && marks[pos+1]) {
            pos++;
            curblockstart = pos;
            while (pos < readupto && marks[pos]) pos++;
            lastblockend = pos;
            qrymove = pos - curblockstart;
            srcmove = correspondance[pos - 1] - correspondance[curblockstart];
            slope = ((double)srcmove) / ((double)qrymove);
            i = curblockstart - 1;
            expected_delta = slope;
            extension_tol = 1 + wiggle_room;
            while (i > 0 && !marks[i] && expected_delta<(double)(correspondance[curblockstart]) && (size_t)expected_delta<maxgap) {
                srcmove = correspondance[curblockstart] - correspondance[i];
                if (correspondance[i] < correspondance[curblockstart] && (double)srcmove > expected_delta/extension_tol && (double)srcmove < expected_delta*extension_tol) {
                    marks[i] = 1;
                }
                i--;
                expected_delta = slope*((double)(curblockstart - i));
            }

            i = lastblockend;
            expected_delta = slope;
            while (i < readupto && !marks[i] && expected_delta<(double)(sourcelen-correspondance[lastblockend-1]) && (size_t)expected_delta<maxgap) {
                srcmove = correspondance[i] - correspondance[lastblockend-1];
                if (correspondance[i] > correspondance[lastblockend - 1] && (double)srcmove > expected_delta/extension_tol && (double)srcmove < expected_delta*extension_tol) {
                    marks[i] = 1;
                }
                i++;
                expected_delta = slope*((double)(i - lastblockend+1));
            }
        }
        pos++;
    }
    pos = 0;
    while (pos < readupto) {
        if (!marks[pos]) correspondance[pos] = sourcelen;
        if (correspondance[pos] >= sourcelen) marks[pos] = 0;
        pos++;
    }
}

char* suftree_approximatesearch_coarse(char* data, size_t datalen, suftree_t* tree, size_t seedalign, size_t minlen, size_t maxgaperror, int add_unaligned_lengths) {
    char* marks;
    size_t* correspondance;
    size_t* blockstarts;
    size_t* blocklens;
    size_t* blockpotential;
    size_t* blockweight;
    size_t numblocks, numblocks_alloc, maxnumblocks;
    size_t curgap,curblock;
    size_t pos, i, j;
    size_t querymove, refmove;
    if (data && datalen > 0) {
        correspondance = malloc(sizeof(size_t)*datalen);
        marks = suftree_markpresent(data, datalen, tree, seedalign/2, correspondance);
        numblocks = 0;
        numblocks_alloc = 2;
        blockstarts = (size_t*) malloc(sizeof(size_t)*numblocks_alloc);
        blocklens = (size_t*)malloc(sizeof(size_t)*numblocks_alloc);
        curblock = 0;
        if (marks[0]) {
            blockstarts[0] = 0;
            numblocks++;
        }
        /* identify block locations */
        for (pos = 1;pos < datalen;pos++) {
            if (marks[pos]) {
                if (!marks[pos] || correspondance[pos - 1] + 2 < correspondance[pos - 1] || correspondance[pos - 1] > correspondance[pos]) {
                    if (numblocks >= numblocks_alloc - 1) {
                        numblocks_alloc <<= 1;
                        blockstarts = (size_t*)realloc(blockstarts, sizeof(size_t)*numblocks_alloc);
                        blocklens = (size_t*)realloc(blocklens, sizeof(size_t)*numblocks_alloc);
                    }
                    blockstarts[numblocks] = pos;
                    curblock = numblocks;
                    numblocks++;
                    blocklens[curblock] = 0;
                }
                blocklens[curblock]++;
            }
        }
        blockpotential = (size_t*)malloc(sizeof(size_t)*numblocks);
        blockweight = (size_t*)malloc(sizeof(size_t)*numblocks);
        /* calculate immediate block extension potential */
        for (i = 0;i < numblocks;i++) {
            j = i + 1;
            /* find the last block to which the current block could potentially be extended */
            while (j < numblocks && blockstarts[j] - blockstarts[i] < tree->datalen + maxgaperror)
                j++;
            /* find an extension that actually makes sense */
            do {
                j--;
                querymove = blockstarts[j] - blockstarts[i];
                if (correspondance[blockstarts[j]] > correspondance[blockstarts[i]]) {
                    refmove = correspondance[blockstarts[j]] - correspondance[blockstarts[i]];
                }
                else
                    refmove = 0;
                if (querymove < refmove) curgap = refmove - querymove;
                else curgap = querymove - refmove;
            } while (j > i && curgap > maxgaperror);
            blockpotential[i] = blockstarts[j] - blockstarts[i] + blocklens[j];
        }
        /* remove blocks with too little potential */
        maxnumblocks = numblocks;
        numblocks = 0;
        memset(marks, 0, datalen);
        for (i = 0;i < maxnumblocks;i++) {
            if (blockpotential[i] > minlen) {
                memset(marks + blockstarts[i], 1, blockpotential[i]);
                numblocks++;
            }
        }
        free(blockstarts);
        free(blocklens);
        free(blockpotential);

        /*expand aligned regions to match the original sequence's length*/
        if (add_unaligned_lengths) {
            while (pos < datalen) {
                if (!marks[pos] && marks[pos - 1]) {
                    j = tree->datalen - correspondance[pos - 1] - 1;
                    if (pos + j > datalen) j = datalen - pos;
                    memset(marks + pos, 1, j);
                    pos += j;
                }
                else if (marks[pos] && !marks[pos - 1]) {
                    j = correspondance[pos];
                    if (j > pos) j = pos;
                    memset(marks + pos - j, 1, j);
                }
                pos++;
            }
        }
        free(correspondance);

        return marks;
    }
    return NULL;
}
char* suftree_approximatesearch_rough(char* data, size_t datalen, suftree_t* tree, size_t seedalign, size_t minlen, size_t maxgaperror, int add_unaligned_lengths) {
    char* marks;
    size_t* correspondance;
    size_t removed, pos, prevblock2, prevblock1, curblock;
    size_t j;
    size_t prevblock2pos, prevblock1pos, curblockpos, prevblocklen, curblocklen;
    int regionswitch;
    /* timing */
#ifdef _DEBUG
    clock_t ct;
    clock_t dct;
    int nf;
    size_t maxj;
    static ht64_t* ht = NULL;
    if (!ht) ht = ht64_alloc_size(100);
    if (0) {
        size_t** names;
        int64_t* values;
        maxj = ht64_astables(ht, (char***)(&names), NULL, &values);
        for (j = 0;j < maxj;j++) free(names[j]);
        free(names);
        free(values);
    }
#endif
    /* get all alignment instances */
    correspondance = malloc(sizeof(size_t)*datalen);
#ifdef _DEBUG
    ct = clock();
#endif
    marks = suftree_markpresent(data, datalen, tree, seedalign/2+1, correspondance);
#ifdef _DEBUG
    dct = clock() - ct;
    ht64_set(ht, (char*)(&datalen), sizeof(size_t), (int64_t)dct, &nf);
#endif
    bridge_substitutions(marks, correspondance, data, datalen, tree->data, tree->datalen, seedalign, 0.2, maxgaperror);
    bridge_substitutions(marks, correspondance, data, datalen, tree->data, tree->datalen, seedalign, 0.2, maxgaperror);
    bridge_substitutions(marks, correspondance, data, datalen, tree->data, tree->datalen, seedalign, 1.0, seedalign);
    /* recursively extend selection over acceptable non-overlapped regions */
    removed = 1;
    while (0) {
        removed = 0;
        prevblock1pos = 0;
        prevblock2pos = 0;
        curblockpos = 0;
        prevblocklen = 0;
        prevblock1 = datalen;
        prevblock2 = datalen;
        curblock = datalen;
        regionswitch = 0;
        curblocklen = 0;
        for (pos = 0;pos < datalen;pos++) {
            if (regionswitch == 0 && marks[pos]) {
                if (!acceptablegap(correspondance, curblockpos + curblocklen - 1, pos, maxgaperror, minlen)) {
                    prevblock2 = prevblock1;
                    prevblock2pos = prevblock1pos;
                    prevblock1 = curblock;
                    prevblock1pos = curblockpos;
                    curblockpos = pos;
                    if (acceptablegap(correspondance, prevblock2pos + prevblocklen - 1, curblockpos, maxgaperror, minlen)) {
                        for (j = prevblock2pos + 1;j < pos;j++) {
                            if (!marks[j])
                                marks[j] = 1;
                        }
                        removed++;
                        curblockpos = prevblock2pos;
                    }
                    else {
                        prevblocklen = curblocklen;
                    }
                }
                else {
                    memset(marks + curblockpos + curblocklen, 1, pos - (curblockpos + curblocklen));
                    curblocklen = pos - curblockpos;
                    removed++;
                }
                regionswitch = 1;
            }
            else if (regionswitch == 1 && !marks[pos]) {
                curblock = correspondance[pos - 1];
                curblocklen = pos - curblockpos;
                regionswitch = 0;
            }
        }
    }
    /* remove marked regions which are too short */
    regionswitch = 0;
    j = 0;
    for (pos = 0;pos < datalen;pos++) {
        if (marks[pos])j++;
        if (regionswitch == 0 && marks[pos]) {
            j = 1;
            regionswitch = 1;
        }
        else if (regionswitch == 1 && marks[pos]==0) {
            if (j < minlen)
                memset(marks + pos - j, 0, j);
            regionswitch = 0;
        }
    }
    if (j < minlen)
        memset(marks + pos - j, 0, j);
    /*expand aligned regions to match the original sequence's length*/
    if (correspondance) {
        pos = 1;
        while (pos < datalen) {
            if (!marks[pos] && marks[pos - 1]) {
                j = tree->datalen - correspondance[pos - 1] - 1;
                if (pos + j > datalen) j = datalen - pos;
                memset(marks + pos, 1, j);
                pos += j;
            }
            else if (marks[pos] && !marks[pos - 1]) {
                j = correspondance[pos];
                if (j > pos) j = pos;
                memset(marks + pos - j, 1, j);
            }
            pos++;
        }
        free(correspondance);
    }
    return marks;
}

char* suftree_approximatesearch(char* data, size_t datalen, suftree_t* tree, size_t seedalign, size_t minlen, size_t maxgaperror, int add_unaligned_lengths) {
    return suftree_approximatesearch_rough(data, datalen, tree, seedalign, minlen, maxgaperror, add_unaligned_lengths);
}

size_t suftree_datalen(suftree_t* tree) {
    return tree->datalen;
}

char** marks_split(char* data, char* marks, size_t datalen, size_t** output_lens, size_t* output_amount, int tokeep) {
    char** result;
    size_t* tmp_output_lens;
    size_t nseqs;
    size_t nseqsalloc;
    size_t i, i0;
    size_t count;
    int mode;
    int inverse;
    nseqs = 0;
    result = NULL;
    mode = 0;
    nseqsalloc = 1;

    tmp_output_lens = (size_t*)calloc(nseqsalloc, sizeof(size_t));
    result = (char**)calloc(nseqsalloc, sizeof(char*));
    inverse = (tokeep != 0);
    i0 = 0;
    if (tokeep == KEEP_UNMARKED_INTERVALS) {
        while (i0 < datalen && !marks[i0]) i0++;
    }
    count = 0;
    for (i = i0;i < datalen;i++) {
        if ( (marks[i]!=0) != mode) {
            mode = (marks[i]!=0);
            if (mode!=inverse) {
                if (nseqs+1 > nseqsalloc) {
                    nseqsalloc *= 2;
                    tmp_output_lens = (size_t*)realloc(tmp_output_lens, sizeof(size_t)*nseqsalloc);
                    result = (char**)realloc(result, sizeof(char*)*nseqsalloc);
                }
            }
            else {
                tmp_output_lens[nseqs] = count;
                result[nseqs] = (char*)malloc(count);
                memcpy(result[nseqs], data + i - count,count);
                nseqs++;
                count = 0;
            }
        }
        if (mode != inverse)count++;
    }
    if (tokeep != KEEP_UNMARKED_INTERVALS && mode != inverse) {
        tmp_output_lens[nseqs] = count;
        result[nseqs] = (char*)malloc(count);
        memcpy(result[nseqs], data + datalen - count, count);
        nseqs++;
    }
    if (nseqs > 0) {
        result = (char**)realloc(result, sizeof(char*)*nseqs);
        tmp_output_lens = (size_t*)realloc(tmp_output_lens, sizeof(size_t)*nseqs);
    }
    else {
        free(result);
        free(tmp_output_lens);
        result = NULL;
        tmp_output_lens = NULL;
    }
    *output_lens = tmp_output_lens;
    *output_amount = nseqs;
    return result;
}
char* mark_correspondance_points(char* querymarks, size_t* correspondance, size_t querylen, size_t sourcelen) {
    char* result;
    size_t i;
    result = calloc(sourcelen, sizeof(char));

    for (i = 0;i < querylen;i++) {
        if (querymarks[i]) {
            result[correspondance[i]] = querymarks[i];
        }
    }
    return result;
}
char* mark_correspondance_intervals(char* querymarks, size_t* correspondance, size_t querylen, size_t sourcelen) {
    /* note: this assumes that the all values between the edges of an interval should be marked. */
    /* important: values within intervals should be increasing */
    char* result;
    size_t i;
    size_t begin;
    int mode;
    result = calloc(sourcelen, sizeof(char));
    mode = 0;
    for (i = 0;i < querylen;i++) {
        if ((querymarks[i] != 0) != mode) {
            if (!mode) {
                begin = correspondance[i];
                mode = 1;
            }
            else {
                memset(result + begin, 1, correspondance[i-1] - begin + 1);
                mode = 0;
            }
        }
    }
    return result;
}
void remove_disordered_marks(char* marks, size_t* correspondance, size_t datalen) {
    int64_t* part_lens;
    size_t* part_starts;
    size_t nseqs;
    size_t nseqsalloc;
    size_t i,j,curseq;
    size_t prevvalue;
    int64_t count;
    size_t curstart;
    size_t log;
    size_t maxcorrespondance;
    int mode;
    int discarded, disordered;
    nseqs = 0;
    mode = 0;
    nseqsalloc = 1;

    part_lens = (int64_t*)calloc(nseqsalloc, sizeof(int64_t));
    part_starts = (size_t*)calloc(nseqsalloc, sizeof(size_t));
    count = 0;
    maxcorrespondance = 0;
    for (i = 0;i < datalen;i++) {
        if (nseqs + 1 > nseqsalloc) {
            nseqsalloc *= 2;
            part_lens = (int64_t*)realloc(part_lens, sizeof(int64_t)*nseqsalloc);
            part_starts = (size_t*)realloc(part_starts, sizeof(size_t)*nseqsalloc);
        }
        if (correspondance[i] > maxcorrespondance) {
            maxcorrespondance = correspondance[i];
            if (i == datalen - 1)maxcorrespondance=correspondance[i]+1;
        }
        if ((marks[i] !=0 ) != mode) {
            mode = (marks[i] != 0);
            if (mode) {
                curstart = i;
            }
            else {
                part_lens[nseqs] = count;
                part_starts[nseqs] = curstart;
                nseqs++;
                count = 0;
            }
        }
        else if (mode && i > 0 && correspondance[i - 1] + 1 != correspondance[i]) {
            part_lens[nseqs] = count;
            part_starts[nseqs] = curstart;
            nseqs++;
            count = 0;
            curstart = i;
        }
        if (mode)count++;
    }
    if (mode) {
        part_lens[nseqs] = count;
        part_starts[nseqs] = curstart;
        nseqs++;
    }
    vec_sort_byi64(part_starts, sizeof(size_t), part_lens, nseqs);
    discarded = 1;
    for (j = 0;j < nseqs ;j++) {
        curseq = nseqs - j - 1;
        disordered = 0;
        if (marks[part_starts[curseq]]) {
            discarded = 0;
            prevvalue = correspondance[part_starts[curseq]];
            log = prevvalue;
            for (i = 0 ;i < part_starts[curseq];i++) {
                if (marks[part_starts[curseq]-i-1]) {
                    if (correspondance[part_starts[curseq] - i - 1] > log) disordered = 1;
                    log = correspondance[part_starts[curseq] - i - 1];
                    if (log > prevvalue) {
                        marks[part_starts[curseq] - i - 1] = 0;
                        correspondance[part_starts[curseq] - i - 1] = maxcorrespondance;
                        discarded = 1;
                    }
                }
            }
            for (i = part_starts[curseq] + 1;i < datalen;i++) {
                if (marks[i]) {
                    if (correspondance[i] < log) disordered = 1;
                    log = correspondance[i];
                    if (log < prevvalue) {
                        marks[i] = 0;
                        correspondance[i] = maxcorrespondance;
                        discarded = 1;
                    }
                }
            }
        }
        if (!disordered)
            break;
    }
    if (part_lens) free(part_lens);
    if (part_starts) free(part_starts);
}
void expand_marks(char* data, char* marks, size_t* correspondance, size_t datalen, char* reference, size_t reflen, int errorcost, int matchgain, int threshold) {
    size_t i;
    size_t j;
    size_t k;
    size_t refi;
    int errorweight;
    int renumber;
    i = 1;
    while (i < datalen) {
        if (marks[i] != marks[i - 1]) {
            if (marks[i]) {
                errorweight = 0;
                refi = correspondance[i];
                j = 1;
                while (errorweight < threshold && refi >= j && i >= j) {
                    if (reference[refi - j] == data[i - j]) {
                        errorweight -= matchgain;
                        marks[i - j] = 1;
                        correspondance[i - j] = refi - j;
                    }
                    else errorweight += errorcost;
                    if (errorweight < 0) errorweight = 0;
                    j++;
                }
                renumber = 0;
                k = j;
                while (k > 0) {
                    k--;
                    if (marks[i - k]) renumber = 1;
                    if (renumber && !marks[i - k])correspondance[i - k] = correspondance[i - k - 1] + 1;
                }
                i++;

            }
            else {
                errorweight = 0;
                refi = correspondance[i-1]+1;
                j = 1;
                while (errorweight < threshold && refi+j<reflen && i+j < datalen) {
                    if (reference[refi + j] == data[i + j]) {
                        errorweight -= matchgain;
                        marks[i + j] = 1;
                        correspondance[i + j] = refi + j;
                    }
                    else errorweight += errorcost;
                    if (errorweight < 0) errorweight = 0;
                    j++;
                }
                renumber = 0;
                k = j;
                while (k > 0) {
                    k--;
                    if (marks[i + k]) renumber = 1;
                    if (renumber && !marks[i + k])correspondance[i + k] = correspondance[i + k + 1] - 1;
                }
                i = i + j;

            }
        }
        else
            i++;
    }
}
size_t count_marks(char* marks, size_t datalen) {
    size_t result;
    size_t i;
    result = 0;
    for (i = 0;i < datalen;i++) {
        if (marks[i])result++;
    }
    return result;
}

size_t* find_high_deltas(size_t* function, size_t datalen, size_t* ngaps, size_t high_delta_threshold, int expected_change) {
    size_t i;
    size_t _ngaps;
    size_t* result;
    long long hdt;

    hdt = (long long)high_delta_threshold;
    _ngaps = 0;
    for (i = 1;i < datalen;i++) {
        if (function[i] > function[i - 1]) {
            if ((long long)(function[i] - function[i - 1]) > hdt + expected_change) _ngaps++;
        }
        else if ((long long)(function[i-1] - function[i]) > hdt - expected_change) _ngaps++;
    }
    result = malloc(sizeof(size_t)*_ngaps);
    *ngaps = _ngaps;
    _ngaps = 0;
    for (i = 1;i < datalen;i++) {
        if (function[i] > function[i - 1]) {
            if ((long long)(function[i] - function[i - 1]) > hdt + expected_change) {
                result[_ngaps] = i;
                _ngaps++;
            }
        }
        else if ((long long)(function[i - 1] - function[i]) > hdt - expected_change) {
            result[_ngaps] = i;
            _ngaps++;
        }
    }
    return result;
}

#define EDIT_INSERTION  0
#define EDIT_DELETION   1
#define EDIT_REPLACE    2
#define EDIT_IDENTITY   3
#define EDIT_NONE       4

static inline size_t edit_cost(int A, int B, int prevedit) {
    if (A == B) return 0;
    if (A == -1) {
        if (prevedit == EDIT_INSERTION || prevedit == EDIT_NONE) return 1;
        return 3;
    }
    if (B == -1) {
        if (prevedit == EDIT_DELETION || prevedit == EDIT_NONE) return 1;
        return 3;
    }
    return 2;
}

size_t min_edit_cost_normal(size_t* r, size_t* identities, size_t* substitutions, int* prev, size_t i, size_t j, size_t p, char* qry, char* ref) {
    size_t cost_insert, cost_delete, cost_replace;
    size_t fc, id, sub;
    /*    O O     .
    .      \|     .
    .     O-?    */
    cost_replace = r[j] + edit_cost(qry[i], ref[i + j - p], prev[j]);
    cost_insert = r[j - 1] + edit_cost(-1, ref[i + j - p], prev[j-1]);
    cost_delete = r[j + 1] + edit_cost(qry[i], -1, prev[j+1]);
    fc = cost_replace;
    id = identities[j];
    sub = substitutions[j];
    if (qry[i] == ref[i + j - p]) {
        id++;
        prev[j] = EDIT_IDENTITY;
    }
    else {
        sub++;
        prev[j] = EDIT_REPLACE;
    }
    if (cost_insert < fc) {
        fc = cost_insert;
        id = identities[j - 1];
        sub = substitutions[j - 1];
        prev[j] = EDIT_INSERTION;
    }
    if (cost_delete < fc) {
        fc = cost_delete;
        id = identities[j+1];
        sub = substitutions[j+1];
        prev[j] = EDIT_DELETION;
    }
    r[j] = fc;
    identities[j] = id;
    substitutions[j] = sub;
    return fc;
}
size_t min_edit_cost_leftmost(size_t* r, size_t* identities, size_t* substitutions, int* prev, size_t i, size_t jmax, size_t p, char* qry, char* ref) {
    /* note that "leftmost" has precedence over "rightmost", */
    /*   which means that a cell which is both will be considered as "leftmost"  */
    size_t cost_delete, cost_replace, j0;
    size_t fc, id, sub;
    if (i + 1 < p) {
        /*      O     .
        .       |     .
        .       ?    */
        j0 = p - i - 1;
        cost_delete = r[j0 + 1] + edit_cost(qry[i], -1, prev[j0+1]);
        fc = cost_delete;
        id = identities[j0 + 1];
        sub = substitutions[j0 + 1];
        prev[j0] = EDIT_DELETION;
    }
    else if (jmax == 1) {
        /*    O       .
        .      \      .
        .       ?    */
        j0 = 0;
        cost_replace = r[0] + edit_cost(qry[i], ref[i - p],prev[0]);
        fc = cost_replace;
        id = identities[0];
        sub = substitutions[0];
        if (qry[i] == ref[i - p]) {
            id++;
            prev[j0] = EDIT_IDENTITY;
        }
        else {
            sub++;
            prev[j0] = EDIT_REPLACE;
        }
    }
    else {
        /*    O O     .
        .      \|     .
        .       ?    */
        j0 = 0;
        cost_replace = r[0] + edit_cost(qry[i], ref[i - p], prev[0]);
        cost_delete = r[1] + edit_cost(qry[i], -1, prev[1]);
        if (cost_delete < cost_replace) {
            fc = cost_delete;
            id = identities[1];
            sub = substitutions[1];
            prev[j0] = EDIT_DELETION;
        }
        else {
            fc = cost_replace;
            id = identities[0];
            sub = substitutions[0];
            if (qry[i] == ref[i - p]) {
                id++;
                prev[j0] = EDIT_IDENTITY;
            }
            else {
                sub++;
                prev[j0] = EDIT_REPLACE;
            }
        }
    }
    r[j0] = fc;
    identities[j0] = id;
    substitutions[j0] = sub;
    return fc;
}
size_t min_edit_cost_rightmost(size_t* r, size_t* identities, size_t* substitutions, int* prev, size_t i, size_t j, size_t p, char* qry, char* ref) {
    /* note that "leftmost" has precedence over "rightmost", */
    /*   which means that a cell which is both will be considered as "leftmost"  */
    size_t cost_insert, cost_replace;
    size_t fc, id, sub;
    /*    O       .
    .      \      .
    .     O-?    */
    cost_replace = r[j] + edit_cost(qry[i], ref[i + j - p], prev[j]);
    cost_insert = r[j - 1] + edit_cost(-1, ref[i + j - p], prev[j-1]);
    if (cost_insert < cost_replace) {
        fc = cost_insert;
        id = identities[j - 1];
        sub = substitutions[j - 1];
        prev[j] = EDIT_INSERTION;
    }
    else {
        fc = cost_replace;
        id = identities[j];
        sub = substitutions[j];
        if (qry[i] == ref[i + j - p]) {
            id++;
            prev[j] = EDIT_IDENTITY;
        }
        else {
            sub++;
            prev[j] = EDIT_REPLACE;
        }
    }
    r[j] = fc;
    identities[j] = id;
    substitutions[j] = sub;
    return fc;
}

size_t _identities_from_edit_distance__shortquery(char* qry, size_t qrylen, char* ref, size_t reflen, size_t* output_alignment_length) {
    size_t identities;
    size_t t, p, offcenter, i, j, j0, jmax, jmaxlocal;
    size_t* r;
    size_t* previdentities;
    size_t* prevsubstitutions;
    int* prevtype;
    size_t delta;
    size_t finalcost;
    size_t diagmotion;
    int acceptable_t;

    offcenter = reflen - qrylen;

    delta = edit_cost(-1, 0, EDIT_INSERTION);
    t = (offcenter + 1) * delta;
    acceptable_t = 0;
    r = NULL;
    previdentities = NULL;
    prevsubstitutions = NULL;
    prevtype = NULL;
    identities = 0;
    while (!acceptable_t && t<reflen*delta*2) {
        p = (size_t)((t / delta) - offcenter)/2;
        if (p > qrylen) p = qrylen;
        jmax = offcenter + 2 * p + 1;
        r = (size_t*)realloc(r, sizeof(size_t)*jmax);
        previdentities = (size_t*)realloc(previdentities, sizeof(size_t)*jmax);
        prevsubstitutions = (size_t*)realloc(prevsubstitutions, sizeof(size_t)*jmax);
        prevtype = (int*)realloc(prevtype, sizeof(int)*jmax);
        r[p] = 0;
        previdentities[p] = 0;
        prevsubstitutions[p] = 0;
        prevtype[p] = EDIT_NONE;
        /* first row */
        for (j = p + 1; j < jmax; j++) {
            r[j] = r[j - 1] + edit_cost(-1, ref[j - p - 1], prevtype[j-1]);
            previdentities[j] = 0;
            prevsubstitutions[j] = 0;
        }
        /* subsequent rows */
        finalcost = t+1;
        for (i = 0;i < qrylen;i++) {
            min_edit_cost_leftmost(r, previdentities, prevsubstitutions, prevtype, i, jmax, p, qry, ref);
            if (i + 1 < p) j0 = p - i - 1;
            else j0 = 0;
            if (jmax + i -p > reflen)jmaxlocal = reflen - i + p;
            else jmaxlocal = jmax - 1;
            for (j = j0 + 1; j < jmaxlocal; j++) {
                min_edit_cost_normal(r, previdentities, prevsubstitutions, prevtype, i, j, p, qry, ref);
            }
            if (jmaxlocal == jmax - 1 && jmax > 1) {
                min_edit_cost_rightmost(r, previdentities, prevsubstitutions, prevtype, i, j, p, qry, ref);
            }
            else
                j--;
            finalcost = r[j];
            identities = previdentities[j];
            if (output_alignment_length) {
                diagmotion = previdentities[j] + prevsubstitutions[j];
                *output_alignment_length = qrylen + reflen - diagmotion;
            }
        }
        acceptable_t = (finalcost <= t);
        t <<= 1;
    }
    if (r)
        free(r);
    if (previdentities)
        free(previdentities);
    if (prevsubstitutions)
        free(prevsubstitutions);
    if (prevtype)
        free(prevtype);
    return identities;
}

size_t identities_from_edit_distance(char* seq1, size_t len1, char* seq2, size_t len2, size_t* output_alignment_length) {
    /* Based on Ukkonen, 1985, Algorithms for approximate string matching */
    if (len1 > len2) return _identities_from_edit_distance__shortquery(seq2, len2, seq1, len1, output_alignment_length);
    else return _identities_from_edit_distance__shortquery(seq1, len1, seq2, len2, output_alignment_length);
}


double suftree_roughalign(char* data, size_t datalen, suftree_t* tree, size_t minmatch, size_t smallgapsize, int trim_mismatches, size_t* alnlen) {
    char* marks;
    size_t* correspondance;
    size_t* indexjumps;
    char* misaligned_query;
    char* misaligned_reference;
    size_t misaligned_query_len;
    size_t misaligned_reference_len;
    size_t ngaps, extragap;
    size_t matches;
    size_t totallen;
    size_t i, j;
    size_t left, right, rleft, rright, index;
    size_t correct_end;
    int count_fragment;
    correspondance = calloc(datalen, sizeof(size_t));
    marks = suftree_markpresent(data, datalen, tree, minmatch, correspondance);
    remove_disordered_marks(marks, correspondance, datalen);
    expand_marks(data, marks, correspondance, datalen, tree->data, tree->datalen, 1, 1, 1);
    matches = 0;/*count_marks(marks, datalen); totallen = datalen;*/
    bridgegaps_1D(marks, datalen, smallgapsize, 1);
    remove_disordered_marks(marks, correspondance, datalen);

    totallen = 0;
    indexjumps = find_high_deltas(correspondance, datalen, &ngaps, 1, 0);
    correct_end = 0;
    i = datalen - 1;
    while (i>0 && correspondance[i] >= tree->datalen - 1) {
        if (correspondance[i - 1] == tree->datalen - 1) correct_end = i;
        i--;
    }
    if (correct_end != 0) {
        if (correct_end > indexjumps[ngaps - 1]) {
            ngaps++;
            indexjumps = realloc(indexjumps,ngaps*sizeof(size_t));
        }
        indexjumps[ngaps - 1] = correct_end;
    }
    i = 0;
    right = left = rright = rleft = 0;
    if (marks[0]) {
        rright = correspondance[0];
        if (!trim_mismatches)
            totallen += correspondance[0];
    }
    if (trim_mismatches)count_fragment = 0;
    else count_fragment = 1;
    while(i<ngaps){
        index = indexjumps[i];
        if (!marks[index - 1] && marks[index]) {
            right = index;
            rright = correspondance[index];
            misaligned_query = data + left;
            misaligned_query_len = right - left;
            misaligned_reference = tree->data + rleft;
            misaligned_reference_len = rright - rleft;
            if (count_fragment) {
                if (rright > rleft + 1) {
                    matches += identities_from_edit_distance(misaligned_query, misaligned_query_len, misaligned_reference, misaligned_reference_len, &extragap);
                    totallen += extragap;
                }
                else if(right - left > 0){
                    totallen += right-left;
                }
                else {
                    totallen++;
                }
            }
        }
        if (marks[index - 1]) {
            left = index;
            rleft = correspondance[index-1]+1;
            if (left + rright > rleft + right) {
                /* the above comparison is equivalent to l-r > rl-rr, but avoids negative numbers */
                totallen += left - right;
            }
            else {
                /* since left is always >= right, in this case "rleft > rright" is guaranteed */
                totallen += rleft - rright;
            }
            for (j = right;j < left;j++) {
                if (data[j] == tree->data[correspondance[j]]) matches++;
            }
            count_fragment = 1;
            if (marks[index]) {
                /* complete deletion between perfect matches */
                right = index;
                rright = correspondance[index];
                totallen += rright - rleft;
            }
        }
        i++;
    }
    if (marks[datalen - 1]) {
        left = datalen;
        rleft = correspondance[datalen-1]+1;
        if (left + rright > rleft + right) {
            /* the above comparison is equivalent to l-r > rl-rr, but avoids negative numbers */
            totallen += left - right;
        }
        else {
            /* since left is always >= right, in this case "rleft > rright" is guaranteed */
            totallen += rleft - rright;
        }
        for (j = right;j < left;j++) {
            if (data[j] == tree->data[correspondance[j]]) matches++;
        }
        if (rleft < tree->datalen && !trim_mismatches)
            totallen += tree->datalen - rleft;
    }
    else if(!trim_mismatches) {
        if (right > left)left = right;
        rleft = correspondance[left - 1] + 1;
        index = datalen-1;
        right = index;
        rright = tree->datalen;
        misaligned_query = data + left;
        misaligned_query_len = right - left;
        misaligned_reference = tree->data + rleft;
        misaligned_reference_len = rright - rleft;
        if (rright > rleft + 1) {
            matches += identities_from_edit_distance(misaligned_query, misaligned_query_len, misaligned_reference, misaligned_reference_len, &extragap);
            totallen += extragap;
        }
        else {
            totallen += right - left;
        }
    }
    free(indexjumps);
    free(marks);
    free(correspondance);
    if (alnlen)*alnlen = totallen;
    
    return 1-(double)matches/(double)totallen;
}

char* suftree_firstsharedregions(suftree_t* reference, suftree_t* query, size_t minlen) {
    suftree_t* tmp;
    tmp = suftree_intersection(reference, query);
    return NULL;
}