#include "datalist.h"
#include <stdlib.h>

int64_t ___added_values___(int64_t delta) {
    static int64_t av = 0;
    av += delta;
    return av;
};

dlink64_t* dlink_alloc(int64_t value){
    dlink64_t* result;
    result = (dlink64_t*)malloc(sizeof(dlink64_t));
    result->value = value;
    result->next = NULL;
    result->prev = NULL;
    return result;
}
dlink64_t* dlink_insert_after(dlink64_t* at, int64_t value) {
    dlink64_t* result;
    ___added_values___(1);
    result = (dlink64_t*)malloc(sizeof(dlink64_t));
    result->value = value;
    result->prev = at;

    if (at) {
        result->next = at->next;
        if (at->next)
            at->next->prev = result;
        at->next = result;
    }
    else
        result->next = NULL;
    return result;
}
dlink64_t* dlink_insert_before(dlink64_t* at, int64_t value) {
    dlink64_t* result;
    ___added_values___(1);
    result = (dlink64_t*)malloc(sizeof(dlink64_t));
    result->value = value;
    result->next = at;
    if (at) {
        result->prev = at->prev;
        if (at->prev)
            at->prev->next = result;
        at->prev = result;
    }
    else
        result->prev = NULL;
    return result;
}
dlink64_t* dlink_append(dlink64_t* list, int64_t value) {
    dlink64_t* result;
    ___added_values___(1);
    result = dlink_alloc(value);
    if (!list) return result;
    while (list->next)
    {
        list = list->next;
    }
    list->next = result;
    result->prev = list;
    return result;
}
dlink64_t* dlink_remove(dlink64_t* target) {
    ___added_values___(-1);
    dlink64_t* result;
    if (target->prev)
        target->prev->next = target->next;
    if (target->next)
        target->next->prev = target->prev;
    if (target->next) result = target->next;
    else if (target->prev) result = target->prev;
    else result = NULL;
    free(target);
    return result;
}
void dlink_freeall(dlink64_t* target) {
    while (target)
        target = dlink_remove(target);
}


hlink64_t* hlink_alloc(int64_t value) {
    hlink64_t* result;
    result = (hlink64_t*)malloc(sizeof(hlink64_t));
    result->value = value;
    result->lchild = NULL;
    result->rchild = NULL;
    result->parent = NULL;
    return result;
}
hlink64_t* hlink_insert(hlink64_t* root, int64_t value) {
    if (!root) {

    }
}
hlink64_t* hlink_remove(hlink64_t* target) {

}
void hlink_freeall(hlink64_t* target) {

}