#include <stdio.h>
#include <math.h> // for NAN
#include <stdlib.h> // for MALLOC
#include <string.h> // for strncpy

// NOTES
// im going to use struct tree_node and a linked list a->b->c->d, 
// e.g. bc merges into a single node a->bc->d where bc points to b and c
// then abc merges into a single node abc->d where abc points to a and bc (which points to b and c)
// then abcd merges into a single node abcd, which is returned.
// i can use this for
// -> MSA!!
// -> gene sequences??

#define IS_ROOT -1

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
    double height;
    int size;
    struct tree_node *next;
    struct tree_node *left;
    struct tree_node *right;
};

struct cell_loc {
    int row;
    int col;
};

struct tree_node *create_list(int elem_no);
struct tree_node *phylo_tree_func(double *dist_matrix, struct tree_node *head, int elem_no);
struct cell_loc shortest_dist_func(double *dist_matrix, int elem_no);

struct tree_node *merge_nodes(struct tree_node *head, int elem_a, int elem_b, int node_num, double dist);
void update_dist_matrix(double *dist_matrix, int elem_a, int elem_b, int size_a, int size_b, int elem_no);
void print_matrix(double *dist_matrix, int elem_no);
void print_list(struct tree_node *head);
void print_tree(struct tree_node *current, int depth, int elem_no);
char *build_full_newick(struct tree_node *head, int elem_no);
void build_newick_recursive(struct tree_node *node, char *newick, int *pos, double parent_height);

// NOTES
// return in newick format
// free memory!!
// error handling

int main(int argc, char *argv[]) {
    char *newick;

    if (argc < 2) {
        fprintf(stderr, "ERROR: Not enough arguments");
        return 1;
    }

    int elem_no = atoi(argv[1]);
    if (argc != 2 + (elem_no * elem_no)) {
        fprintf(stderr, "ERROR: Not enough values in matrix");
        return 1;
    }

    double *dist_matrix = malloc(elem_no * elem_no * sizeof(double));
    
    // dist_matrix is 1d so I can pass values to main easier
    int arg_index = 2;
    for (int i = 0; i < elem_no; i++) {
        for (int j = 0; j < elem_no; j++) {
            dist_matrix[i * elem_no + j] = atof(argv[arg_index]);
            arg_index++;
        }
    }

    struct tree_node *head = create_list(elem_no);

    fprintf(stderr, "\n------ INITIAL MATRIX ------\n");
    print_matrix(dist_matrix, elem_no);

    fprintf(stderr, "\n------ START UPGMA ------\n");

    head = phylo_tree_func(dist_matrix, head, elem_no);

    fprintf(stderr, "\n------ FINAL TREE ------\n");
    print_tree(head, 0, elem_no);

    fprintf(stderr, "\n------ FINAL NEWICK ------\n");
    newick = build_full_newick(head, elem_no);
    printf("%s", newick);

    return 0;
}

struct tree_node *create_list(int elem_no) {
    struct tree_node *head = NULL;
    struct tree_node *current = head;
    struct tree_node *new;

    for (int i = 0; i < elem_no; i++) {
        new = malloc(sizeof(struct tree_node));
        new->dist_mat_id = i;
        new->node_id = i;
        new->height = 0;
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
    }
    current->next = NULL;

    return head;
}

struct tree_node *phylo_tree_func(double *dist_matrix, struct tree_node *head, int elem_no) {
    struct cell_loc shortest_dist_cell;
    struct tree_node *current;
    int size_a;
    int size_b;
    int node_num = elem_no;
    int finished = FALSE;
    int shortest_dist_row;
    int shortest_dist_col;
    
    fprintf(stderr, "\nINITIAL TREE\n");
    print_list(head);
    print_matrix(dist_matrix, elem_no);

    while (finished == FALSE) {
        shortest_dist_cell = shortest_dist_func(dist_matrix, elem_no);
        shortest_dist_row = shortest_dist_cell.row;
        shortest_dist_col = shortest_dist_cell.col;
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

        head = merge_nodes(head, shortest_dist_row, shortest_dist_col, node_num, dist_matrix[shortest_dist_row * elem_no + shortest_dist_col]);
        node_num++;
        update_dist_matrix(dist_matrix, shortest_dist_row, shortest_dist_col, size_a, size_b, elem_no);
        if (head->next == NULL) {
            finished = TRUE;
        }
        
        fprintf(stderr, "\nMERGE NO. %d\n", node_num - elem_no);
        print_list(head);
        print_matrix(dist_matrix, elem_no);
    }

    return head;

}

struct cell_loc shortest_dist_func(double *dist_matrix, int elem_no) {
    struct cell_loc cell_loc;
    double cell;
    double min = NAN;

    for (int row = 0; row < elem_no; row++) {
        for (int col = 0; col < elem_no; col++) {
            if (row != col && row < elem_no && col < elem_no ) {
                cell = dist_matrix[row * elem_no + col];
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


struct tree_node *merge_nodes(struct tree_node *head, int elem_a, int elem_b, int node_num, double dist) {
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
    new_node->height = dist/2;
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

void update_dist_matrix(double *dist_matrix, int elem_a, int elem_b, int size_a, int size_b, int elem_no) {
    double new_dist;

    for (int i = 0; i < elem_no; i++) {
        if (i != elem_a && i != elem_b && !isnan(dist_matrix[elem_a * elem_no + i]) && !isnan(dist_matrix[elem_b * elem_no + i])) {
            new_dist = (dist_matrix[i * elem_no + elem_a] * size_a + dist_matrix[elem_b * elem_no + i] * size_b) / (size_a + size_b);
            dist_matrix[i * elem_no + elem_a] = new_dist;
            dist_matrix[elem_a * elem_no + i] = new_dist;
        }
    }

    // replace second element values with NAN
    for (int j = 0; j < elem_no; j++) {
        dist_matrix[j * elem_no + elem_b] = NAN;
        dist_matrix[elem_b * elem_no + j] = NAN;
    }
}

void print_matrix(double *dist_matrix, int elem_no) {
    for (int i = 0; i < elem_no * 6 + 2; i++) {
        fprintf(stderr, "-");
    }
    for (int i = 0; i < elem_no; i++) {
        fprintf(stderr, "\n||");
        for (int j = 0; j < elem_no; j++) {
            if (!isnan(dist_matrix[i * elem_no + j])) {
                fprintf(stderr, "%.2f||", dist_matrix[i * elem_no + j]);
            } else {
                fprintf(stderr, "    ||");
            }
        }
    }
    fprintf(stderr, "\n");
    for (int i = 0; i < elem_no * 6 + 2; i++) {
        fprintf(stderr, "-");
    }
    fprintf(stderr, "\n");
}

void print_list(struct tree_node *head) {
    struct tree_node *current = head;
    int list_pos = 0;
    int left_index = -1;
    int right_index = -1;
    int next_index = -1;
    int current_index = -1;

    fprintf(stderr, "\n");
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

        fprintf(stderr, "(NODE %d = %d)\n", list_pos, current->node_id);
        fprintf(stderr, "     Left = %d\n", left_index);
        fprintf(stderr, "     Right = %d\n", right_index);
        fprintf(stderr, "     Next = %d\n", next_index);
        
        left_index = -1;
        right_index = -1;
        next_index = -1;

        list_pos++;

        current = current->next;
    }
    fprintf(stderr, "\n");
}

void print_tree(struct tree_node *current, int depth, int elem_no) {
    if (current != NULL) {

        for (int i = 0; i < depth; i++) {
            fprintf(stderr, "|    ");
        }

        if (current->node_id < elem_no) {
            fprintf(stderr, "'--- ");
            fprintf(stderr, "%d\n", current->node_id);
        } else {
            fprintf(stderr, "|--- ");
            fprintf(stderr, "o\n"); // intermediate node - not an organism
        }

        print_tree(current->left, depth + 1, elem_no);
        print_tree(current->right, depth + 1, elem_no);
    }
    return;
}

char *build_full_newick(struct tree_node *head, int elem_no) {
    char *newick = malloc(elem_no * 20); // arbitrary number - ensures adequately large malloc
    int pos = 0;

    build_newick_recursive(head, newick, &pos, IS_ROOT);
    newick[pos] = ';';
    pos++;
    newick[pos] = '\0';
    
    return newick;
}

void build_newick_recursive(struct tree_node *node, char *newick, int *pos, double parent_height) {
    double branch_len;
    int digits_added;

    if (node == NULL) {
        return;
    }

    if (node->left == NULL && node->right == NULL) {
        branch_len = parent_height - node->height;
        digits_added = sprintf(newick + *pos, "%d:%.2f", node->node_id, branch_len);
        *pos += digits_added;
        return;
    }

    newick[*pos] = '(';
    (*pos)++;

    build_newick_recursive(node->left, newick, pos, node->height);
    newick[*pos] = ',';
    (*pos)++;
    build_newick_recursive(node->right, newick, pos, node->height);

    newick[*pos] = ')';
    (*pos)++;

    if (parent_height != IS_ROOT) {
        branch_len = parent_height - node->height;
        digits_added = sprintf(newick + *pos, ":%.2f", branch_len);
        *pos += digits_added;
    }
}