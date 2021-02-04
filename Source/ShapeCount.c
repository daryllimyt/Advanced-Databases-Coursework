#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <assert.h>

typedef void* SortMergeJoinDatabase;
SortMergeJoinDatabase SortMergeJoinAllocateDatabase(unsigned long sizeNumberOfEdgesInTheEnd);
void SortMergeJoinInsertEdge(SortMergeJoinDatabase database, int fromNodeID, int toNodeID, int edgeLabel);
int SortMergeJoinRunQuery(SortMergeJoinDatabase database, int edgeLabel1, int edgeLabel2, int edgeLabel3);
void SortMergeJoinDeleteEdge(SortMergeJoinDatabase database, int fromNodeID, int toNodeID, int edgeLabel);
void SortMergeJoinDeleteDatabase(SortMergeJoinDatabase database);

typedef void* HashjoinDatabase;
HashjoinDatabase HashjoinAllocateDatabase(unsigned long sizeNumberOfEdgesInTheEnd);
void HashjoinInsertEdge(HashjoinDatabase database, int fromNodeID, int toNodeID, int edgeLabel);
int HashjoinRunQuery(HashjoinDatabase database, int edgeLabel1, int edgeLabel2, int edgeLabel3);
void HashjoinDeleteEdge(HashjoinDatabase database, int fromNodeID, int toNodeID, int edgeLabel);
void HashjoinDeleteDatabase(HashjoinDatabase database);

typedef void* CompetitionDatabase;
CompetitionDatabase CompetitionAllocateDatabase(unsigned long sizeNumberOfEdgesInTheEnd);
void CompetitionInsertEdge(CompetitionDatabase database, int fromNodeID, int toNodeID, int edgeLabel);
int CompetitionRunQuery(CompetitionDatabase database, int edgeLabel1, int edgeLabel2, int edgeLabel3);
void CompetitionDeleteEdge(CompetitionDatabase database, int fromNodeID, int toNodeID, int edgeLabel);
void CompetitionDeleteDatabase(CompetitionDatabase database);

typedef void* T;
typedef struct {
    int fromNodeID;
    int toNodeID;
    int edgeLabel;
} CustomTuple;
CustomTuple* customTupleInit(int fromNodeID, int toNodeID, int edgeLabel);

typedef struct {
    T* items;
    unsigned long capacity;
    unsigned long size;
} CustomVector;
CustomVector* customVectorInit(unsigned long capacity);
void customVectorAdd(CustomVector *v, T item); 
void customVectorFree(CustomVector *v);
void customVectorDelete(CustomVector *v, int index);

void swap(T* a, T* b);
int countFinalTable(CustomVector* table);
void freeTables(CustomVector* el1, CustomVector* el2, CustomVector* el3, CustomVector* j1, CustomVector* j2, int** hashTable, unsigned long htSize);
int partition(CustomVector* v, int lo, int hi, int attributeIndex);
void quickSort(CustomVector* v, int lo, int hi, int attributeIndex);
CustomVector* merge(CustomVector* lhs, CustomVector* rhs, int lhsAttributeIndex, int rhsAttributeIndex);
CustomVector* sortMergeVectors(CustomVector* lhs, CustomVector* rhs, int lhsAttributeIndex, int rhsAttributeIndex);
CustomVector* hashJoinVectors(CustomVector* lhs, CustomVector* rhs, int lhsAttributeIndex, int rhsAttributeIndex, int** hashTable, unsigned long htSize);
void clearHashTable(int** hashTable, unsigned long htSize);

// Helper functions and custom data structure functions

CustomTuple* customTupleInit(int fromNodeID, int toNodeID, int edgeLabel){
    CustomTuple* v = (CustomTuple*)malloc(sizeof(CustomTuple));
    assert(v != NULL);
    v->fromNodeID = fromNodeID;
    v->toNodeID = toNodeID;
    v->edgeLabel = edgeLabel;
    return v;
}

CustomVector* customVectorInit(unsigned long capacity){
    CustomVector* v = (CustomVector*)malloc(sizeof(CustomVector));
    assert(v != NULL);
    v->items = (T*)malloc(sizeof(T) * capacity);
    v->capacity = capacity;
    v->size = 0;
    return v;
}

void customVectorAdd(CustomVector *v, T item){
    assert(v->size < v->capacity);
    v->items[v->size++] = item;
}

void customVectorFree(CustomVector *v){
    for (int i = 0; i < v->size; i++) free(v->items[i]);
    free(v->items);
    free(v);
}
void customVectorDelete(CustomVector* v, int index){
    if(index < 0 || index >= (int)v->capacity) return;
    for (int i = index + 1; i < v->size; i++) swap(&v->items[i-1], &v->items[i]);
    free((v->items[--v->size])); // free the last element
    v->capacity--;
}

void swap(T* a, T* b) { 
    T t = *a; 
    *a = *b; 
    *b = t; 
} 

int countFinalTable(CustomVector* table){
    int count = 0;
    int* row;
    for(int i = 0; i < table->size; i++){
        row = (int*)table->items[i];
        if(row[3] == row[0]) {
            count++; //el3.to = el1.from
        }
    }
    return count;
}

void freeTables(CustomVector* el1, CustomVector* el2, CustomVector* el3, CustomVector* j1, CustomVector* j2, int** hashTable, unsigned long htSize){
    if (el1 != NULL) customVectorFree(el1);
    if (el2 != NULL) customVectorFree(el2);
    if (el3 != NULL) customVectorFree(el3);
    if (j1 != NULL) customVectorFree(j1);
    if (j2 != NULL) customVectorFree(j2);
    if (hashTable != NULL && htSize != 0) {
        for(int i = 0; i < htSize; i++){
            if(hashTable[i] != NULL) free(hashTable[i]);
        }
        free(hashTable);
    }
}

// Sort Merge Database functions

CustomVector* sortMergeVectors(CustomVector* lhs, CustomVector* rhs, int lhsAttributeIndex, int rhsAttributeIndex) {
    quickSort(lhs, 0, lhs->size-1, lhsAttributeIndex);
    quickSort(rhs, 0, rhs->size-1, rhsAttributeIndex);
    return merge(lhs,rhs,lhsAttributeIndex,rhsAttributeIndex);
}

SortMergeJoinDatabase SortMergeJoinAllocateDatabase(unsigned long sizeNumberOfEdgesInTheEnd){
    return customVectorInit(sizeNumberOfEdgesInTheEnd);
}

void SortMergeJoinInsertEdge(SortMergeJoinDatabase database, int fromNodeID, int toNodeID, int edgeLabel) {
    customVectorAdd((CustomVector*)database, customTupleInit(fromNodeID, toNodeID, edgeLabel));
}

int SortMergeJoinRunQuery(SortMergeJoinDatabase database, int edgeLabel1, int edgeLabel2, int edgeLabel3) {
    CustomVector* v = (CustomVector*)database; 
    if (v->size < 3) return 0; // Not enough edges to form a triangle
    CustomTuple* t;
    int* row;
    int label1 = 0, label2 = 0, label3 = 0;
    for(int i = 0; i < v->size; i++) {
        t = (CustomTuple*)v->items[i];
        if (t->edgeLabel == edgeLabel1) label1++;
        if (t->edgeLabel == edgeLabel2) label2++;
        if (t->edgeLabel == edgeLabel3) label3++;
    }
    CustomVector* el1 = customVectorInit(label1); // joins 1t = 2f
    CustomVector* el2 = customVectorInit(label2); // joins 2t = 3f
    CustomVector* el3 = customVectorInit(label3); // joins 3t = 1f
    for(int i = 0; i < v->size; i++) {
        // Convert tuple to array
        t = (CustomTuple*)v->items[i];
        // Need to malloc separately to sort each edge label table
        if (t->edgeLabel == edgeLabel1) {
            row = (int*)malloc(sizeof(int)*2);
            row[0] = t->fromNodeID;
            row[1] = t->toNodeID;
            customVectorAdd(el1, row);
        } 
        if (t->edgeLabel == edgeLabel2) {
            row = (int*)malloc(sizeof(int)*2);
            row[0] = t->fromNodeID;
            row[1] = t->toNodeID;
            customVectorAdd(el2, row);
        } 
        if (t->edgeLabel == edgeLabel3) {
            row = (int*)malloc(sizeof(int)*2);
            row[0] = t->fromNodeID;
            row[1] = t->toNodeID;
            customVectorAdd(el3, row);
        } 
    }

    if(el1->size == 0 || el2->size == 0 || el3->size == 0) {
        freeTables(el1, el2, el3, NULL, NULL, NULL, 0);
        return 0;
    }

    CustomVector* j1 = sortMergeVectors(el1, el2, 1, 0);
    if (j1->size == 0) {
        freeTables(el1, el2, el3, j1, NULL, NULL, 0);
        return 0;
    }
    CustomVector* j2 = sortMergeVectors(j1, el3, 2, 0);
 
    int res = countFinalTable(j2);
    freeTables(el1, el2, el3, j1, j2, NULL, 0);
    return res;
}

void SortMergeJoinDeleteEdge(SortMergeJoinDatabase database, int fromNodeID, int toNodeID, int edgeLabel) {
    CustomTuple* t;
    CustomVector* db = (CustomVector*)database;
    for(int i = 0; i < db->size; i++){ //for each item in the database, may be duplicated
        t = (CustomTuple*)db->items[i];
        if (t->fromNodeID == fromNodeID && t->toNodeID == toNodeID && t->edgeLabel == edgeLabel) { //check if everything in tuple is same
            customVectorDelete((CustomVector*)database, i); 
        }
    }
}

void SortMergeJoinDeleteDatabase(SortMergeJoinDatabase database) {
    customVectorFree((CustomVector*)database);
}

int partition(CustomVector* v, int lo, int hi, int attributeIndex){ 
    int* row = (int*)v->items[hi];    // get last row
    int pivot = row[attributeIndex]; // get pivot value from the comparator attribute
    int i = lo-1;  // Index of smaller element 
    
    for (int j = lo; j < hi; j++) { 
        row = (int*)v->items[j]; 
        if (row[attributeIndex] < pivot) swap(&v->items[++i], &v->items[j]);  // Move smaller values to i
    } 
    swap(&v->items[i+1], &v->items[hi]); // Swap pivot row to correct position
    return i+1; 
} 

void quickSort(CustomVector* v, int lo, int hi, int attributeIndex){ 
    if (lo < hi) { 
        int pivotIndex = partition(v, lo, hi, attributeIndex); // Returns correct position of v[pivorIndex]
        quickSort(v, lo, pivotIndex-1, attributeIndex); 
        quickSort(v, pivotIndex+1, hi, attributeIndex); 
    } 
} 

CustomVector* merge(CustomVector* lhs, CustomVector* rhs, int lhsAttributeIndex, int rhsAttributeIndex) {
    CustomVector* res = customVectorInit(lhs->size*rhs->size);
    int i = 0, j = 0, mark = -1; // mark is to help with duplicates
    int *lhsRow, *rhsRow, *row;
    while (i < lhs->size && j < rhs->size) {
        // If mismatch, move pointer
        if (mark < 0) {
            while (i < lhs->size-1 && ((int*)lhs->items[i])[lhsAttributeIndex] < ((int*)rhs->items[j])[rhsAttributeIndex]) i++;
            while (j < rhs->size-1 && ((int*)lhs->items[i])[lhsAttributeIndex] > ((int*)rhs->items[j])[rhsAttributeIndex]) j++;
            mark = j; // To mark first occurrence of duplicate item
        }
        lhsRow = (int*)lhs->items[i];
        rhsRow = (int*)rhs->items[j];

        if (rhsRow[rhsAttributeIndex] == lhsRow[lhsAttributeIndex]) {
            row = (int*)malloc(sizeof(int)*4);
            for (int i = 0; i <= lhsAttributeIndex; i++) row[i] = lhsRow[i];
            row[lhsAttributeIndex+1] = rhsRow[1];
            customVectorAdd(res, row);
            j++;
            if (j >= rhs->size) {
                j = mark;
                i++;
            }
        } else {
            j = mark; // rewind to first occurrence
            i++;
            mark = -1;
        }
    }
    return res;
}

// Hash Join Database functions

HashjoinDatabase HashjoinAllocateDatabase(unsigned long sizeNumberOfEdgesInTheEnd){
    return customVectorInit(sizeNumberOfEdgesInTheEnd);
}

void HashjoinInsertEdge(HashjoinDatabase database, int fromNodeID, int toNodeID, int edgeLabel){
    customVectorAdd((CustomVector*)database, customTupleInit(fromNodeID, toNodeID, edgeLabel)); 
}

int HashjoinRunQuery(HashjoinDatabase database, int edgeLabel1, int edgeLabel2, int edgeLabel3){
    CustomVector* v = (CustomVector*)database;
    if (v->size < 3) return 0; // Not enough edges to form a triangle
    // unsigned long htSize = v->size * 2;
    CustomTuple* t;
    int* row;
    int label1 = 0, label2 = 0, label3 = 0;
    for(int i = 0; i < v->size; i++) {
        t = (CustomTuple*)v->items[i];
        if (t->edgeLabel == edgeLabel1) label1++;
        if (t->edgeLabel == edgeLabel2) label2++;
        if (t->edgeLabel == edgeLabel3) label3++;
    }
    CustomVector* el1 = customVectorInit(label1); // joins 1t = 2f
    CustomVector* el2 = customVectorInit(label2); // joins 2t = 3f
    CustomVector* el3 = customVectorInit(label3); // joins 3t = 1f
    for(int i = 0; i < v->size; i++) {
        // Convert tuple to array
        t = (CustomTuple*)v->items[i];
        
        // Need to malloc separately to sort each edge label table
        if (t->edgeLabel == edgeLabel1) {
            row = (int*)malloc(sizeof(int)*2);
            row[0] = t->fromNodeID;
            row[1] = t->toNodeID;
            customVectorAdd(el1, row);
        } 
        if (t->edgeLabel == edgeLabel2) {
            row = (int*)malloc(sizeof(int)*2);
            row[0] = t->fromNodeID;
            row[1] = t->toNodeID;
            customVectorAdd(el2, row);
        } 
        if (t->edgeLabel == edgeLabel3) {
            row = (int*)malloc(sizeof(int)*2);
            row[0] = t->fromNodeID;
            row[1] = t->toNodeID;
            customVectorAdd(el3, row);
        } 
    }

    if(el1->size == 0 || el2->size == 0 || el3->size == 0) {
        freeTables(el1,el2,el3,NULL,NULL,NULL,0);
        return 0;
    }
    unsigned long htSize;
    if (label1 >= label2 && label1 >= label3) {
        htSize = 2 * label1;
    } else if (label2 >= label1 && label2 >= label3) {
        htSize = 2 * label2;
    } else {
        htSize = 2 * label3;
    }
    int** hashTable = (int**)calloc(htSize,sizeof(int*));
    // join on [el1.from, el1.to] [el2.from, el2.to] el1.to = el2.from
    CustomVector* j1 = hashJoinVectors(el1, el2, 1, 0, hashTable, htSize); 
    if(j1->size == 0){
        freeTables(el1,el2,el3,j1,NULL,hashTable,htSize);
        return 0;
    }
    clearHashTable(hashTable, htSize);
    // join on [el1.from, el1.to=el2.from, el2.to] [el3.from, el3.to] el2.to = el3.from
    CustomVector* j2 = hashJoinVectors(j1, el3, 2, 0, hashTable, htSize); 
    int value = countFinalTable(j2);
    freeTables(el1,el2,el3,j1,j2,hashTable,htSize);
    return value;
}

void HashjoinDeleteEdge(HashjoinDatabase database, int fromNodeID, int toNodeID, int edgeLabel){
    CustomTuple* t;
    CustomVector* db = (CustomVector*)database;
    for(int i = 0; i < db->size; i++){ //for each item in the database, may be duplicated
        t = (CustomTuple*)db->items[i];
        //check if everything in tuple matches edge
        if (t->fromNodeID == fromNodeID && t->toNodeID == toNodeID && t->edgeLabel == edgeLabel) { 
            customVectorDelete((CustomVector*)database, i); 
        }
    }
}

void HashjoinDeleteDatabase(HashjoinDatabase database){
    customVectorFree((CustomVector*)database);
}

CustomVector* hashJoinVectors(CustomVector* lhs, CustomVector* rhs, int lhsAttributeIndex, int rhsAttributeIndex, int** hashTable, unsigned long htSize){
    // always want RHS to be smaller -> hashtable built using RHS elements, then probe using LHS
    // if the condition below is met, then order is swapped. lhs rightmost == rhs leftmost -> lhs leftmost == rhs rightmost
    // the table that is being joined on the 0th index has row len == 2
    // for the other table, the element is being joined on the rightmost element.
    if (lhs->size < rhs->size) {
        return hashJoinVectors(rhs, lhs, rhsAttributeIndex, lhsAttributeIndex, hashTable, htSize);
    }
    CustomVector* res = customVectorInit(lhs->size*rhs->size);
    
    // Build HT using RHS
    int rhsRowLen = (rhsAttributeIndex == 0) ? 2 : rhsAttributeIndex + 1;
    int lhsRowLen = (lhsAttributeIndex == 0) ? 2 : lhsAttributeIndex + 1;
    for(int i = 0; i < rhs->size; i++){
        int* buildInput = (int*)rhs->items[i];
        int hashValue = buildInput[rhsAttributeIndex] % htSize; 
        while(hashTable[hashValue] != NULL){
            hashValue = ++hashValue % htSize;
        }
        int* htRow = (int*)malloc(sizeof(int)*rhsRowLen); // hash array row size == rhs width
        for (int j = 0; j < rhsRowLen; j++) {
            htRow[j] = buildInput[j];
        }
        hashTable[hashValue] = htRow; 
    }

    // Probe HT using LHS
    int newRowLen = rhsRowLen + lhsRowLen - 1; // Either 3 or 4
    for(int i = 0; i < lhs->size; i++){
        int* probeInput = (int*)lhs->items[i];
        int hashValue = probeInput[lhsAttributeIndex] % htSize; 
        int* htRow = (int*)hashTable[hashValue];
        while(htRow != NULL){
            if(htRow[rhsAttributeIndex] == probeInput[lhsAttributeIndex]){
                int* newRow = (int*)malloc(sizeof(int)*newRowLen);
                if (rhsAttributeIndex == 0 && lhsAttributeIndex == 1 && newRowLen == 3) {
                    newRow[0] = probeInput[0]; //el1.from
                    newRow[1] = probeInput[1]; //el1.to = el2.from
                    newRow[2] = htRow[1]; //el2.to
                } else if(rhsAttributeIndex == 1 && lhsAttributeIndex == 0 && newRowLen == 3){
                    newRow[0] = htRow[0]; //el1.from
                    newRow[1] = htRow[1]; //el1.to = el2.from
                    newRow[2] = probeInput[1]; //el2.to
                } else if(rhsAttributeIndex == 0 && lhsAttributeIndex == 2 && newRowLen == 4){
                    newRow[0] = probeInput[0]; //el1.from 
                    newRow[1] = probeInput[1]; //el1.to = el2.from 
                    newRow[2] = htRow[0]; //el3.from
                    newRow[3] = htRow[1]; //el3.to
                } else if(rhsAttributeIndex == 2 && lhsAttributeIndex == 0 && newRowLen == 4){
                    newRow[0] = htRow[0]; //el1.from
                    newRow[1] = htRow[1]; //el1.to = el2.from 
                    newRow[2] = probeInput[0]; //el3.from
                    newRow[3] = probeInput[1]; //el3.to
                }
                customVectorAdd(res,newRow);
            }
            hashValue = ++hashValue % htSize;
            htRow = hashTable[hashValue];
        }
    }
    return res;
}

void clearHashTable(int** hashTable, unsigned long htSize){
    for(int i = 0; i < htSize; i++){
        if(hashTable[i] != NULL) {
            free(hashTable[i]);
            hashTable[i] = NULL;
        }
    }
}

// Competition Database functions

CompetitionDatabase CompetitionAllocateDatabase(unsigned long sizeNumberOfEdgesInTheEnd){
    return customVectorInit(sizeNumberOfEdgesInTheEnd);
}

void CompetitionInsertEdge(CompetitionDatabase database, int fromNodeID, int toNodeID, int edgeLabel){
    customVectorAdd((CustomVector*)database, customTupleInit(fromNodeID, toNodeID, edgeLabel)); 
}

int CompetitionRunQuery(CompetitionDatabase database, int edgeLabel1, int edgeLabel2, int edgeLabel3){
    return SortMergeJoinRunQuery(database, edgeLabel1, edgeLabel2, edgeLabel3);
}

void CompetitionDeleteEdge(CompetitionDatabase database, int fromNodeID, int toNodeID, int edgeLabel){
    CustomTuple* t;
    CustomVector* db = (CustomVector*)database;
    for(int i = 0; i < db->size; i++){ //for each item in the database, may be duplicated
        t = (CustomTuple*)db->items[i];
        //check if everything in tuple matches edge
        if (t->fromNodeID == fromNodeID && t->toNodeID == toNodeID && t->edgeLabel == edgeLabel) { 
            customVectorDelete((CustomVector*)database, i); 
        }
    }
}

void CompetitionDeleteDatabase(CompetitionDatabase database){
    customVectorFree((CustomVector*)database);
}
