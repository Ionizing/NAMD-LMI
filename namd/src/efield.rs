use pest_derive::Parser;
use shared::Result;


mod fnparse {
    use once_cell::sync::Lazy;
    use pest_derive::Parser;
    use pest::{
        Parser,
        iterators::Pairs,
        pratt_parser::{
            Assoc,
            Op,
            PrattParser,
        },
    };

    #[derive(Parser)]
    #[grammar_inline = r#"
    value = @{ int ~ ("." ~ ASCII_DIGIT*)? ~ (^"e" ~ int)? }
        int = { ("+" | "-")? ~ ASCII_DIGIT+ }
    infix = _{ add | sub | mul | div | pow }
        add = { "+" }
        sub = { "-" }
        mul = { "*" }
        div = { "/" }
        pow = { "^" }
    prefix = _{ neg | exp }
        neg = { "-" }
        exp = { "e^" }
    function = _{ sin | cos | tan }
        sin = { "sin" }
        cos = { "cos" }
        tan = { "tan" }
    variable = @{ "t" }
    primary = _{
        variable |
        value |
        function ~ "(" ~ expr ~ ")" |
        "(" ~ expr ~ ")"
    }
    expr = {
        prefix? ~ primary ~ (infix ~ prefix? ~ primary)*
    }
    program = _{ SOI ~ expr ~ EOI }
    WHITESPACE = _{ " " | "\t" | "\n" | "\r" }
    "#]
    struct Function;

    static PRATT_PARSER: Lazy<PrattParser<Rule>> = Lazy::new(|| {
        use Rule::*;
        use Assoc::*;

        PrattParser::new()
            .op(Op::infix(add, Left) | Op::infix(sub, Right))
            .op(Op::infix(mul, Left) | Op::infix(div, Right))
            .op(Op::prefix(sin) | Op::prefix(cos) | Op::prefix(tan))
            .op(Op::prefix(neg))
            .op(Op::infix(pow, Right))
            .op(Op::prefix(exp))
    });


    #[derive(Debug)]
    pub enum Operation {
        // binary op
        Add,
        Sub,
        Mul,
        Div,
        Pow,

        // unary op
        Neg,
        Exp,

        Sin,
        Cos,
        Tan,
    }


    #[derive(Debug)]
    pub enum Expr {
        Variable,
        Value(f64),
        UnaryOp {
            op: Operation,
            rhs: Box<Expr>,
        },
        BinOp {
            lhs: Box<Expr>,
            op: Operation,
            rhs: Box<Expr>,
        },
    }


    pub fn parse_expr(pairs: Pairs<Rule>) -> Expr {
        PRATT_PARSER
            .map_primary(|primary| match primary.as_rule() {
                Rule::variable => Expr::Variable,
                Rule::value    => Expr::Value(primary.as_str().parse::<f64>().unwrap()),
                Rule::expr     => parse_expr(primary.into_inner()),
                _ => unreachable!("Expr::parse Expected primary expression, found {:?}", primary),
            })
            .map_prefix(|op, rhs| {
                let op = match op.as_rule() {
                    Rule::neg => Operation::Neg,
                    Rule::exp => Operation::Exp,
                    Rule::sin => Operation::Sin,
                    Rule::cos => Operation::Cos,
                    Rule::tan => Operation::Tan,
                    _ => unreachable!("Expr::parse Expected prefix operator, found {:?}", op),
                };
                Expr::UnaryOp {
                    op,
                    rhs: Box::new(rhs),
                }
            })
            .map_infix(|lhs, op, rhs| {
                let op = match op.as_rule() {
                    Rule::add => Operation::Add,
                    Rule::sub => Operation::Sub,
                    Rule::mul => Operation::Mul,
                    Rule::div => Operation::Div,
                    Rule::pow => Operation::Pow,
                    _ => unreachable!("Expr::parse Expected prefix operator, found {:?}", op),
                };
                Expr::BinOp {
                    lhs: Box::new(lhs),
                    op,
                    rhs: Box::new(rhs),
                }
            })
            .parse(pairs)
    }


    impl Expr {
        fn eval(&self, var: f64) -> f64 {
            match self {
                Expr::Variable => var,
                Expr::Value(x) => *x,
                Expr::UnaryOp{op, rhs} => {
                    match op {
                        Operation::Neg => - rhs.eval(var),
                        Operation::Exp => rhs.eval(var).exp(),
                        Operation::Sin => rhs.eval(var).sin(),
                        Operation::Cos => rhs.eval(var).cos(),
                        Operation::Tan => rhs.eval(var).tan(),
                        _ => unreachable!()
                    }
                }
                Expr::BinOp{lhs, op, rhs} => {
                    match op {
                        Operation::Add => lhs.eval(var) + rhs.eval(var),
                        Operation::Sub => lhs.eval(var) - rhs.eval(var),
                        Operation::Mul => lhs.eval(var) * rhs.eval(var),
                        Operation::Div => lhs.eval(var) / rhs.eval(var),
                        Operation::Pow => lhs.eval(var).powf(rhs.eval(var)),
                        _ => unreachable!()
                    }
                }
            }
        }
    }


    pub fn str2fn(i: &str) -> impl Fn(f64) -> f64 {
        let mut pairs = Function::parse(Rule::program, i).unwrap();
        let expr = parse_expr(pairs.next().unwrap().into_inner());
        move |x| expr.eval(x)
    }
}


/// Example of input external-field data:
/// ```ignore
/// # comments
/// lcycle = false # or true
///
/// # Vector3 data
/// efield = Vector3 {
/// #    Ex   Ey   Ez
///     0.0  0.0  0.0
///     1.0  2.0  3.0
///     ...
/// }
///
/// # or three Function (Ex, Ey, Ez)
/// efield = Function3 {
///     sin(t) * e^(- 0.001 * t^2 + 500)
///     cos(t) * e^(- 0.001 * t^2 + 500)
///     0
/// }
/// ```
#[derive(Parser)]
#[grammar_inline = r##"
float = @{ int ~ ("." ~ ASCII_DIGIT*)? ~ (^"e" ~ int)? }
    int = _{ ("+" | "-")? ~ ASCII_DIGIT+ }

vector3 = {
    vector3tag ~ "{" ~
        (float ~ float ~ float) ~
        ("," ~ float ~ float ~ float)* ~ (",")? ~
    "}"
}
    vector3tag = @{ ^"vector3" }

function3 = {
    functiontag ~ "{" ~
        functionbody ~ ";" ~
        functionbody ~ ";" ~
        functionbody ~ (";")? ~
    "}"
}
    functiontag  = @{ ^"function3" }
    functionbody = @{ not_semicolon_or_brace+ }
    not_semicolon_or_brace = {
        !(
            ";" |
            "{" |
            "}"
        )
        ~ ANY
    }

efield     = _{ SOI ~ ( vector3 | function3 ) ~ EOI }
WHITESPACE = _{ " " | "\n" | "\t" | "\r" }
COMMENT    = _{ "#" ~ (!NEWLINE ~ ANY)* ~ NEWLINE }
"##]
pub enum Efield {
    Vector3(Vec<[f64; 3]>),
    Function3([Box<dyn Fn(f64) -> f64>; 3]),
}


impl Efield {
    fn string_to_vec3(s: &str) -> Vec<[f64; 3]> {
        let efield_vec = s
            .split(&[',', ' ', '\n', '\t', '\r'])
            .filter(|x| !x.is_empty())
            .map(|v| v.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();
        
        assert_eq!(efield_vec.len() % 3, 0,
            "Length of EFIELD data should be in multiples of 3, got {}.", efield_vec.len());

        efield_vec.chunks_exact(3)
            .map(|v| [v[0], v[1], v[2]])
            .collect::<Vec<_>>()
    }

    fn string_to_function(s: &str) -> Box<dyn Fn(f64) -> f64> {
        Box::new(fnparse::str2fn(s))
    }

    pub fn from_str(s: &str) -> Result<Self> {
        use pest::Parser;
        use pest::iterators::Pair;

        fn parse_internal(pair: Pair<Rule>) -> Efield {
            match pair.as_rule() {
                Rule::vector3   => {
                    let token = pair.as_str()
                        .split(&['{', '}'])
                        .nth(1)
                        .unwrap();
                    Efield::Vector3(Efield::string_to_vec3(token))
                },
                Rule::function3 => {
                    let token = pair.as_str()
                        .split(&['{', '}'])
                        .nth(1)
                        .unwrap();
                    let tokens = token.split(";")
                        .map(|x| x.trim())
                        .collect::<Vec<&str>>();
                    Efield::Function3([
                        Efield::string_to_function(tokens[0]),
                        Efield::string_to_function(tokens[1]),
                        Efield::string_to_function(tokens[2]),
                    ])
                },
                _ => panic!("Unexpected token. Only function3 and vector3 are acceptable.")
            }
        }

        let efield = Efield::parse(Rule::efield, s).unwrap().next().unwrap();
        Ok(parse_internal(efield))
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_string_to_vec3() {
        let s = r#"
        1 2 3,
        3 2 1 "#;
        let vec3 = Efield::string_to_vec3(s);
        assert_eq!(vec3, &[[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]);
    }

    #[test]
    fn test_string_to_function() {
        let s = "1 + sin(t) * cos(t) + e^(-t)";
        let f = Efield::string_to_function(s);
        assert_eq!(2.0, f(0.0));
        assert!( (1.8225281546 - f(1.0)) < 1e-8 );
    }

    #[test]
    fn test_parse_efield_function3() {
        let s = r#"
        Function3 {
            e^(-0.001 * (t - 500)^2) * sin(t);
            e^(-0.001 * (t - 500)^2) * cos(t);
            0
        }"#;
        match Efield::from_str(s).unwrap() {
            Efield::Function3(f) => {
                assert!((f[0](498.0) - 0.9943582286).abs() < 1E-8);
                assert!((f[1](498.0) + 0.05730294897).abs() < 1E-8);
                assert!((f[2](498.0) - 0.0).abs() < 1E-8);
            },
            _ => panic!("Parse failed.")
        }
    }


    #[test]
    fn test_parse_efield_vector3() {
        let s = r#"
        Vector3 {
            1 2 3,
            3 2 1,
        }"#;
        match Efield::from_str(s).unwrap() {
            Efield::Vector3(vec3) => {
                assert_eq!(vec3, &[[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]);
            },
            _ => panic!("Parse failed.")
        }
    }
}
