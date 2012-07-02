#ifndef IMAGESTACK_PARSER_H
#define IMAGESTACK_PARSER_H

#include "Image.h"
#include "Statistics.h"
#include "header.h"

#ifndef roundf
#define roundf(x) floorf(x+0.5);
#endif

class Expression {
    // the AST nodes are below here
    /* Grammar:

    IfThenElse -> Condition | Condition ? Condition : Condition
    Condition  -> Sum > Sum | Sum < Sum |
    Sum >= Sum | Sum <= Sum |
    Sum == Sum | Sum != Sum | Sum
    Sum     -> Sum + Product | Sum - Product | Product
    Product -> Product * Factor | Product / Factor |
    Product % Factor | Factor
    Factor  -> Term ^ Term | Term
    Term    -> Funct0 ( ) | Funct1 ( IfThenElse ) | Funct2 ( IfThenElse , IfThenElse ) | - Term | Var | ( IfThenElse ) | Float | Sample | Const | Uniform
    Funct0  -> mean | sum | max | min | stddev | var | skew | kurtosis
    Funct1  -> sin | cos | tan | log | abs | mean | sum | max | min | stddev | variance | skew | kurtosis
    Funct2  -> covariance
    Var     -> x | y | t | c | val
    Uniform -> width | height | frames | channels
    Const   -> e | pi
    Sample  -> [IfThenElse, IfThenElse, IfThenElse] | [IfThenElse, IfThenElse] | [IfThenElse]

    */
public:

    struct State {
        State(Image im_) : x(0), y(0), t(0), c(0), im(im_), stats(im_) {}
        int x, y, t, c;
        Image im;
        Stats stats;
    };

    struct Node {
        Node() {};
        virtual ~Node() {};
        virtual float eval(State &state) = 0;
    };

    struct Unary : public Node {
        Unary(Node *arg_) : arg(arg_) {}
        ~Unary() {delete arg;}
        Node *arg;
    };

    struct Binary : public Node {
        Binary(Node *left_, Node *right_) : left(left_), right(right_) {}
        ~Binary() {delete left; delete right;}
        Node *left, *right;
    };

    struct Ternary : public Node {
        Ternary(Node *left_, Node *middle_, Node *right_) : left(left_), middle(middle_), right(right_) {}
        ~Ternary() {delete left; delete middle; delete right;}
        Node *left, *middle, *right;
    };

    struct Negation : public Unary {
        Negation(Node *a) : Unary(a) {}
        float eval(State &state) {return -arg->eval(state);}
    };

    struct IfThenElse : public Ternary {
        IfThenElse(Node *l, Node *m, Node *r) : Ternary(l, m, r) {}
        float eval(State &state) {return left->eval(state) ? middle->eval(state) : right->eval(state);}
    };

    struct LTE : public Binary {
        LTE(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) <= right->eval(state) ? 1 : 0;}
    };

    struct GTE : public Binary {
        GTE(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) >= right->eval(state) ? 1 : 0;}
    };

    struct LT : public Binary {
        LT(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) < right->eval(state) ? 1 : 0;}
    };

    struct GT : public Binary {
        GT(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) > right->eval(state) ? 1 : 0;}
    };

    struct EQ : public Binary {
        EQ(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) == right->eval(state) ? 1 : 0;}
    };

    struct NEQ : public Binary {
        NEQ(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) != right->eval(state) ? 1 : 0;}
    };

    struct Plus : public Binary {
        Plus(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) + right->eval(state);}
    };

    struct Minus : public Binary {
        Minus(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) - right->eval(state);}
    };

    struct Mod : public Binary {
        Mod(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return fmod(left->eval(state), right->eval(state));}
    };

    struct Times : public Binary {
        Times(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) * right->eval(state);}
    };

    struct Divide : public Binary {
        Divide(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return left->eval(state) / right->eval(state);}
    };

    struct Power : public Binary {
        Power(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return powf(left->eval(state), right->eval(state));}
    };

    struct Funct_sin : public Unary {
        Funct_sin(Node *a) : Unary(a) {}
        float eval(State &state) {return sinf(arg->eval(state));}
    };

    struct Funct_cos : public Unary {
        Funct_cos(Node *a) : Unary(a) {}
        float eval(State &state) {return cosf(arg->eval(state));}
    };

    struct Funct_tan : public Unary {
        Funct_tan(Node *a) : Unary(a) {}
        float eval(State &state) {return tanf(arg->eval(state));}
    };

    struct Funct_atan : public Unary {
        Funct_atan(Node *a) : Unary(a) {}
        float eval(State &state) {return atanf(arg->eval(state));}
    };

    struct Funct_asin : public Unary {
        Funct_asin(Node *a) : Unary(a) {}
        float eval(State &state) {return asinf(arg->eval(state));}
    };

    struct Funct_acos : public Unary {
        Funct_acos(Node *a) : Unary(a) {}
        float eval(State &state) {return acosf(arg->eval(state));}
    };

    struct Funct_atan2 : public Binary {
        Funct_atan2(Node *l, Node *r) : Binary(l, r) {}
        float eval(State &state) {return atan2f(left->eval(state), right->eval(state));}
    };

    struct Funct_abs : public Unary {
        Funct_abs(Node *a) : Unary(a) {}
        float eval(State &state) {return fabsf(arg->eval(state));}
    };

    struct Funct_floor : public Unary {
        Funct_floor(Node *a) : Unary(a) {}
        float eval(State &state) {return floorf(arg->eval(state));}
    };

    struct Funct_ceil : public Unary {
        Funct_ceil(Node *a) : Unary(a) {}
        float eval(State &state) {return ceilf(arg->eval(state));}
    };

    struct Funct_round : public Unary {
        Funct_round(Node *a) : Unary(a) {}
        float eval(State &state) {return roundf(arg->eval(state));}
    };

    struct Funct_log : public Unary {
        Funct_log(Node *a) : Unary(a) {}
        float eval(State &state) {return logf(arg->eval(state));}
    };

    struct Funct_exp : public Unary {
        Funct_exp(Node *a) : Unary(a) {}
        float eval(State &state) {return expf(arg->eval(state));}
    };

    struct Funct_mean0 : public Node {
        float eval(State &state) {return state.stats.mean();}
    };

    struct Funct_mean1 : public Unary {
        Funct_mean1(Node *a) : Unary(a) {}
        float eval(State &state) {return state.stats.mean((int)(arg->eval(state) + 0.5));}
    };

    struct Funct_sum0 : public Node {
        float eval(State &state) {return state.stats.sum();}
    };

    struct Funct_sum1 : public Unary {
        Funct_sum1(Node *a) : Unary(a) {}
        float eval(State &state) {return state.stats.sum((int)(arg->eval(state) + 0.5));}
    };

    struct Funct_max0 : public Node {
        float eval(State &state) {return state.stats.maximum();}
    };

    struct Funct_max1 : public Unary {
        Funct_max1(Node *a) : Unary(a) {}
        float eval(State &state) {return state.stats.maximum((int)(arg->eval(state) + 0.5));}
    };

    struct Funct_min0 : public Node {
        float eval(State &state) {return state.stats.minimum();}
    };

    struct Funct_min1 : public Unary {
        Funct_min1(Node *a) : Unary(a) {}
        float eval(State &state) {return state.stats.minimum((int)(arg->eval(state) + 0.5));}
    };

    struct Funct_variance0 : public Node {
        float eval(State &state) {return state.stats.variance();}
    };

    struct Funct_variance1 : public Unary {
        Funct_variance1(Node *a) : Unary(a) {}
        float eval(State &state) {return state.stats.variance((int)(arg->eval(state) + 0.5));}
    };

    struct Funct_stddev0 : public Node {
        float eval(State &state) {return sqrtf(state.stats.variance());}
    };

    struct Funct_stddev1 : public Unary {
        Funct_stddev1(Node *a) : Unary(a) {}
        float eval(State &state) {return sqrtf(state.stats.variance((int)(arg->eval(state) + 0.5)));}
    };

    struct Funct_skew0 : public Node {
        float eval(State &state) {return state.stats.skew();}
    };

    struct Funct_skew1 : public Unary {
        Funct_skew1(Node *a) : Unary(a) {}
        float eval(State &state) {return state.stats.skew((int)(arg->eval(state) + 0.5));}
    };

    struct Funct_kurtosis0 : public Node {
        float eval(State &state) {return state.stats.kurtosis();}
    };

    struct Funct_kurtosis1 : public Unary {
        Funct_kurtosis1(Node *a) : Unary(a) {}
        float eval(State &state) {return state.stats.kurtosis((int)(arg->eval(state) + 0.5));}
    };

    struct Funct_covariance : public Binary {
        Funct_covariance(Node *left_, Node *right_) : Binary(left_, right_) {}
        float eval(State &state) {return state.stats.covariance((int)(left->eval(state) + 0.5), (int)(right->eval(state) + 0.5));}
    };

    struct SampleHere : public Unary {
        SampleHere(Node *a) : Unary(a) {}
        float eval(State &state) {
            int c = (int)(arg->eval(state) + 0.5);
            return state.im(state.x, state.y, state.t, c);
        }
    };

    struct Sample2D : public Binary {
        Sample2D(Node *left_, Node *right_) : Binary(left_, right_) {
        }

        float eval(State &state) {
            if (sample.size() != (size_t)state.im.channels) sample.resize(state.im.channels);
            state.im.sample2D(left->eval(state), right->eval(state), state.t, sample);
            return sample[state.c];
        }

        vector<float> sample;
    };

    struct Sample3D : public Ternary {
        Sample3D(Node *left_, Node *middle_, Node *right_) : Ternary(left_, middle_, right_) {
        }

        float eval(State &state) {
            if (sample.size() != (size_t)state.im.channels) sample.resize(state.im.channels);
            state.im.sample3D(left->eval(state),
                              middle->eval(state),
                              right->eval(state), sample);
            return sample[state.c];
        }

        vector<float> sample;
    };

    struct Var_x : public Node {
        float eval(State &state) {return state.x;}
    };

    struct Var_y : public Node {
        float eval(State &state) {return state.y;}
    };

    struct Var_t : public Node {
        float eval(State &state) {return state.t;}
    };

    struct Var_c : public Node {
        float eval(State &state) {return state.c;}
    };

    struct Var_val : public Node {
        float eval(State &state) {return state.im(state.x, state.y, state.t, state.c);}
    };

    struct Uniform_width : public Node {
        float eval(State &state) {return state.im.width;}
    };

    struct Uniform_height : public Node {
        float eval(State &state) {return state.im.height;}
    };

    struct Uniform_frames : public Node {
        float eval(State &state) {return state.im.frames;}
    };

    struct Uniform_channels : public Node {
        float eval(State &state) {return state.im.channels;}
    };

    struct Float : public Node {
        Float(float value_) : value(value_) {}
        float eval(State &state) {return value;}
        float value;
    };

    // all the parsing stuff is below here
private:
    void skipWhitespace();
    bool match(string prefix);
    bool consume(string prefix);

    // IfThenElse -> Condition (? Condition : Condition)?
    Node *parseIfThenElse();

    // Condition -> Sum ((>|<|>=|<=|==|!=) Sum)?
    Node *parseCondition();

    // Sum     -> Product ((+|-) Product)*
    Node *parseSum();

    // Product -> Factor ((*|/|%) Factor)*
    Node *parseProduct();

    // Factor  -> Term ^ Term | Term
    Node *parseFactor();

    // Term    -> Funct ( ) | Funct ( IfThenElse , IfThenElse ) | Funct ( IfThenElse ) | - Term | Var | ( IfThenElse ) | Float | Sample
    Node *parseTerm();

    Node *root;
    string source;
    size_t sourceIndex;
    bool varyingAllowed;

public:
    Expression(string source_, bool varyingAllowed = true);

    ~Expression();

    float eval(State &state);

    static void help();

};

#include "footer.h"
#endif
