//use std::fmt::Display;

use pest_derive::Parser;
use shared::Result;

pub use fnparse::Expr;

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


    #[derive(Debug, Clone)]
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


    #[derive(Debug, Clone)]
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
        pub fn eval(&self, var: f64) -> f64 {
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

        // TODO: better display
        pub fn to_str(&self) -> String {
            match self {
                Expr::Variable => "t".to_string(),
                Expr::Value(x) => x.to_string(),
                Expr::UnaryOp { op, rhs } => {
                    match op {
                        Operation::Neg => format!("- {}", rhs.to_str()),
                        Operation::Exp => format!("e^({})", rhs.to_str()),
                        Operation::Sin => format!("sin({})", rhs.to_str()),
                        Operation::Cos => format!("cos({})", rhs.to_str()),
                        Operation::Tan => format!("tan({})", rhs.to_str()),
                        _ => unreachable!(),
                    }
                }
                Expr::BinOp { lhs, op, rhs } => {
                    match op {
                        Operation::Add => format!("({}) + ({})", lhs.to_str(), rhs.to_str()),
                        Operation::Sub => format!("({}) - ({})", lhs.to_str(), rhs.to_str()),
                        Operation::Mul => format!("({}) * ({})", lhs.to_str(), rhs.to_str()),
                        Operation::Div => format!("({}) / ({})", lhs.to_str(), rhs.to_str()),
                        Operation::Pow => format!("({})^({})", lhs.to_str(), rhs.to_str()),
                        _ => unreachable!(),
                    }
                }
            }
        }
    }


    pub fn str2expr(i: &str) -> Expr {
        let mut pairs = Function::parse(Rule::program, i).unwrap();
        parse_expr(pairs.next().unwrap().into_inner())
    }


    pub fn str2fn(i: &str) -> impl Fn(f64) -> f64 {
        let expr = str2expr(i);
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

vector3 = { vector3tag ~ "{" ~ vector3body ~ "}" }
    vector3tag = @{ ^"vector3" }
    vector3body = {
        (float ~ float ~ float) ~
        ("," ~ float ~ float ~ float)* ~ (",")?
    }

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
#[derive(Clone)]
pub enum Efield {
    Vector3(Vec<[f64; 3]>),
    Function3([Expr; 3]),
}


impl Efield {
    #[allow(dead_code)]
    fn string_to_function(s: &str) -> Box<dyn Fn(f64) -> f64> {
        Box::new(fnparse::str2fn(s))
    }

    pub fn from_str(s: &str) -> Result<Self> {
        use pest::Parser;
        use pest::iterators::Pair;

        fn parse_internal(pair: Pair<Rule>) -> Efield {
            match pair.as_rule() {
                Rule::vector3   => {
                    let vec = pair.clone().into_inner()
                        .skip(1)
                        .next().unwrap()
                        .into_inner()
                        .filter_map(|tok| {
                            if tok.as_rule() == Rule::float {
                                tok.as_str().parse::<f64>().ok()
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<f64>>()
                        .chunks_exact(3)
                        .map(|v| [v[0], v[1], v[2]])
                        .collect();

                    Efield::Vector3(vec)
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
                        fnparse::str2expr(tokens[0]),
                        fnparse::str2expr(tokens[1]),
                        fnparse::str2expr(tokens[2]),
                    ])
                },
                _ => panic!("Unexpected token. Only function3 and vector3 are acceptable.")
            }
        }

        let efield = Efield::parse(Rule::efield, s).unwrap().next().unwrap();
        Ok(parse_internal(efield))
    }


    /// Substitute the 't' with actual time, and get the real-time electric field.
    ///
    /// - For Function3, the result is exact;
    /// - For Vector3, the result is linearly interpolated if the time index is not integer.
    pub fn eval(&self, t: f64, dt: f64) -> [f64; 3] {
        use Efield::*;

        match self {
            Function3([f1, f2, f3]) => {
                [
                    f1.eval(t),
                    f2.eval(t),
                    f3.eval(t),
                ]
            },
            Vector3(vec) => {
                let t = t / dt;

                if t == (t as u32) as f64 {
                    let t = t as usize;
                    return vec[t];
                }

                let t0 = t.floor() as usize;
                let t1 = t0 + 1;
                let tmt0 = t - t0 as f64;
                [
                    (vec[t1][0] - vec[t0][0]) * tmt0 + vec[t0][0],
                    (vec[t1][1] - vec[t0][1]) * tmt0 + vec[t0][1],
                    (vec[t1][2] - vec[t0][2]) * tmt0 + vec[t0][2],
                ]
            },
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;

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

        let efield = Efield::from_str(s).unwrap();
        let evaluated = efield.eval(498.0, 1.0);

        assert!((evaluated[0] - 0.9943582286).abs() < 1E-8);
        assert!((evaluated[1] + 0.05730294897).abs() < 1E-8);
        assert!((evaluated[2] - 0.0).abs() < 1E-8);
    }

    #[test]
    fn test_parse_efield_vector3() {
        let s = r#"
        Vector3 {
            1 2 3,
            3 2 1,
        }"#;

        let efield = Efield::from_str(s).unwrap();
        let evaluated0 = efield.eval(0.0, 1.0);
        let evaluated1 = efield.eval(0.5, 1.0);
        let evaluated2 = efield.eval(1.0, 1.0);

        assert_eq!(evaluated0, [1.0, 2.0, 3.0]);
        assert_eq!(evaluated1, [2.0, 2.0, 2.0]);
        assert_eq!(evaluated2, [3.0, 2.0, 1.0]);
    }

    #[test]
    fn test_parse_efield_vector3_with_comment() {
        let s = r#"
        Vector3 {
        #   x y z  time
            1 2 3, # 0fs
            3 2 1, # 1fs
        }"#;
        match Efield::from_str(s).unwrap() {
            Efield::Vector3(vec3) => {
                assert_eq!(vec3, &[[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]);
            },
            _ => panic!("Parse failed.")
        }
    }
}
