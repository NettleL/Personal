#include <stdio.h>
#include <math.h> // for NAN
#include <stdlib.h> // for MALLOC
#include <string.h> // for strncpy

#define MAX_NAME_LEN 20
#define MAX_ELEM 6

// NOTES
// im going to use struct tree_node and a linked list a->b->c->d, 
// e.g. bc merges into a single node a->bc->d where bc points to b and c
// then abc merges into a single node abc->d where abc points to a and bc (which points to b and c)
// then abcd merges into a single node abcd, which is returned.

enum Active {
    IS_ACTIVE,
    IS_INACTIVE
};

enum Truth {
    TRUE,
    FALSE
};

struct tree_node {
    int dist_mat_id;
    int node_id;
    double branch_len;
    int size;
    struct tree_node *next;
    struct tree_node *left;
    struct tree_node *right;
};

struct cell_loc {
    int row;
    int col;
};

struct dict {
    int index;
    char name[MAX_NAME_LEN];
    struct dict *next;
};

struct dict *create_dict(char *elem_names[], int elem_no);
struct tree_node *create_list(int elem_no);
struct tree_node *phylo_tree_func(double dist_matrix[MAX_ELEM][MAX_ELEM], struct tree_node *head, int elem_no);
struct cell_loc shortest_dist_func(double dist_matrix[MAX_ELEM][MAX_ELEM], int elem_no);

struct tree_node *merge_nodes(struct tree_node *head, int elem_a, int elem_b, int node_num);
void update_dist_matrix(double dist_matrix[MAX_ELEM][MAX_ELEM], int elem_a, int elem_b, int size_a, int size_b);
void print_matrix(double dist_matrix[MAX_ELEM][MAX_ELEM]);
void print_list(struct tree_node *head);
void print_tree(struct tree_node *current, struct dict *dict_head, int depth, int elem_no);
char *get_name(struct dict *head, int index);

int main(void) {
    double dist_matrix[MAX_ELEM][MAX_ELEM] = {
        {0, 1, 2, 3, 5},
        {1, 0, 1, 2, 1},
        {2, 1, 0, 4, 1},
        {3, 2, 4, 0, 1},
        {5, 1, 1, 1, 0}
    };
    char *elem_names[] = {"kitty", "doggie", "snek", "floofball", "ducc"};
    int elem_no = 5;

    struct dict *dict_head = create_dict(elem_names, elem_no);
    struct tree_node *head = create_list(elem_no);

    printf("INIT MATRIX\n");
    print_matrix(dist_matrix);

    printf("START\n");

    head = phylo_tree_func(dist_matrix, head, elem_no);

    printf("\n============== FINAL TREE ==============\n");
    print_tree(head, dict_head, 0, elem_no);

    return 0;
}

struct dict *create_dict(char *elem_names[], int elem_no) {
    struct dict *head = NULL;
    struct dict *current = head;
    struct dict *new;

    for (int i = 0; i < elem_no; i++) {
        new = malloc(sizeof(struct dict));

        strncpy(new->name, elem_names[i], sizeof(new->name) - 1); // is there a better way to do this?? surely
        new->name[sizeof(new->name) - 1] = '\0';
        new->index = i;

        if (head == NULL) {
            head = new;
            current = head;
        } else {
            current->next = new;
            current = current->next;
        }

        printf("(%d: %s) -> ", current->index, current->name);
    }
    printf("\n");
    current->next = NULL;

    return head;
}

struct tree_node *create_list(int elem_no) {
    struct tree_node *head = NULL;
    struct tree_node *current = head;
    struct tree_node *new;

    for (int i = 0; i < elem_no; i++) {
        new = malloc(sizeof(struct tree_node));
        new->dist_mat_id = i;
        new->node_id = i;
        new->branch_len = NAN;
        new->size = 1;
        new->left = NULL;
        new->right = NULL;
        new->next = NULL;

        if (head == NULL) {
            head = new;
            current = head;
        } else {
            current->next = new;
            current = current->next;
        }
        
        printf("%d -> ", current->node_id);
    }
    printf("\n");
    current->next = NULL;

    return head;
}

struct tree_node *phylo_tree_func(double dist_matrix[MAX_ELEM][MAX_ELEM], struct tree_node *head, int elem_no) {
    struct cell_loc shortest_dist_cell;
    struct tree_node *current;
    int size_a;
    int size_b;
    int node_num = elem_no;
    int finished = FALSE;
    
    print_list(head);
    print_matrix(dist_matrix);

    while (finished == FALSE) {
        shortest_dist_cell = shortest_dist_func(dist_matrix, elem_no);

        current = head;
        size_a = 1;
        size_b = 1;

        while (current != NULL) {
            if (current->dist_mat_id == shortest_dist_cell.row) {
                size_a = current->size;
            }
            if (current->dist_mat_id == shortest_dist_cell.col) {
                size_b = current->size;
            }
            current = current->next;
        }

        head = merge_nodes(head, shortest_dist_cell.row, shortest_dist_cell.col, node_num);
        node_num++;
        update_dist_matrix(dist_matrix, shortest_dist_cell.row, shortest_dist_cell.col, size_a, size_b);
        if (head->next == NULL) {
            finished = TRUE;
        }
        
        print_list(head);
        print_matrix(dist_matrix);
    }

    return head;

}

struct cell_loc shortest_dist_func(double dist_matrix[MAX_ELEM][MAX_ELEM], int elem_no) {
    struct cell_loc cell_loc;
    double cell;
    double min = NAN;

    for (int row = 0; row < MAX_ELEM; row++) {
        for (int col = 0; col < MAX_ELEM; col++) {
            if (row != col && row < elem_no && col < elem_no ) {
                cell = dist_matrix[row][col];
                if (!isnan(cell) && (isnan(min) || cell < min)) {
                    cell_loc.row = row;
                    cell_loc.col = col;
                    min = cell;
                }
            }
        }
    }
    return cell_loc;
}


struct tree_node *merge_nodes(struct tree_node *head, int elem_a, int elem_b, int node_num) {
    struct tree_node *elem_i = NULL;
    struct tree_node *elem_j = NULL;

    struct tree_node *current = head;
    struct tree_node *next = NULL;
    struct tree_node *prev = NULL;
    struct tree_node *tail;
    
    struct tree_node *new_node = malloc(sizeof(struct tree_node));

    while (current != NULL) {
        if (current->dist_mat_id == elem_a) {
            elem_i = current;
        }
        if (current->dist_mat_id == elem_b) {
            elem_j = current;
        }
        current = current->next;
    }

    current = head;
    while (current != NULL) {
        next = current->next;
        if (current == elem_i || current == elem_j) {
            if (prev == NULL) {
                head = current->next;
            } else {
                prev->next = current->next;
            }
        } else {
            prev = current;
        }
        current = next;
    }

    elem_i->next = NULL;
    elem_j->next = NULL;
    
    new_node->dist_mat_id = elem_a;
    new_node->node_id = node_num;
    new_node->size = elem_i->size + elem_j->size;
    new_node->left = elem_i;
    new_node->right = elem_j;
    new_node->next = NULL;
    
    tail = head;
    if (head == NULL) {
        head = new_node;
    } else {
        while (tail->next != NULL) {
            tail = tail->next;
        }
        tail->next = new_node;
    }

    return head;
}

void update_dist_matrix(double dist_matrix[MAX_ELEM][MAX_ELEM], int elem_a, int elem_b, int size_a, int size_b) {
    double new_dist;

    for (int i = 0; i < MAX_ELEM; i++) {
        if (i != elem_a && i != elem_b && !isnan(dist_matrix[elem_a][i]) && !isnan(dist_matrix[elem_b][i])) {
            new_dist = (dist_matrix[i][elem_a] * size_a + dist_matrix[elem_b][i] * size_b) / (size_a + size_b);
            dist_matrix[i][elem_a] = new_dist;
            dist_matrix[elem_a][i] = new_dist;
        }
    }

    // replace second element values with NAN
    for (int j = 0; j < MAX_ELEM; j++) {
        dist_matrix[j][elem_b] = NAN;
        dist_matrix[elem_b][j] = NAN;
    }
}

void print_matrix(double dist_matrix[MAX_ELEM][MAX_ELEM]) {
    printf("======================================");
    for (int i = 0; i < MAX_ELEM; i++) {
        printf("\n||");
        for (int j = 0; j < MAX_ELEM; j++) {
            if (!isnan(dist_matrix[i][j])) {
                printf("%.2f||", dist_matrix[i][j]);
            } else {
                printf("    ||");
            }
        }
    }
    printf("\n======================================\n");
}

void print_list(struct tree_node *head) {
    struct tree_node *current = head;
    int list_pos = 0;
    int left_index = -1;
    int right_index = -1;
    int next_index = -1;
    int current_index = -1;

    printf("\n");
    while (current != NULL) {
        if (current->left != NULL) {
            left_index = current->left->node_id;
        } 

        if (current->right != NULL) {
            right_index = current->right->node_id;
        } 

        if (current->next != NULL) {
            next_index = current->next->node_id;
        } 

        printf("======= NODE %d = %d =======\n", list_pos, current->node_id);
        printf("LEFT = %d\n", left_index);
        printf("RIGHT = %d\n", right_index);
        printf("NEXT = %d\n", next_index);
        
        left_index = -1;
        right_index = -1;
        next_index = -1;

        list_pos++;

        current = current->next;
    }
    printf("\n");
}

void print_tree(struct tree_node *current, struct dict *dict_head, int depth, int elem_no) {
    if (current != NULL) {

        for (int i = 0; i < depth; i++) {
            printf("|    ");
        }

        if (current->node_id < elem_no) {
            printf("'--- ");
            printf("%s\n", get_name(dict_head, current->node_id));
        } else {
            printf("|--- ");
            printf("o\n", current->node_id); // intermediate node - not an organism
        }

        print_tree(current->left, dict_head, depth + 1, elem_no);
        print_tree(current->right, dict_head, depth + 1, elem_no);
    }
    return;
}

char *get_name(struct dict *head, int index) {
    struct dict *current = head;
    for (int i = 0; i < index; i++) {
        current = current->next;
    }
    return current->name;
}