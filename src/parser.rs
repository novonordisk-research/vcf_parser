use pom::parser::*;
use std::str::{self, FromStr};
use serde_json::{Value, json};

#[derive(Debug)]
#[derive(serde::Serialize)]
pub enum PropertyVal {
    SimpleValue(Value),
    Group(Vec<Value>),
}

// Define the parser for whitespace
fn space<'a>() -> Parser<'a, u8, ()> {
    one_of(b" \t\r\n").repeat(0..).discard()
}

// Define the parser for identifiers
fn ident<'a>() -> Parser<'a, u8, String> {
    one_of(b"._abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789").repeat(1..).convert(String::from_utf8)
}

// Define the parser for operators
fn operator<'a>() -> Parser<'a, u8, String> {
    (seq(b"==") 
    | seq(b"is") 
    | seq(b"=") 
    | seq(b"<=") 
    | seq(b">=")
    | seq(b"!=")
    | seq(b"<") 
    | seq(b">")
    | seq(b"eq")
    | seq(b"ne")
    | seq(b"gt")
    | seq(b"ge")
    | seq(b"lt")
    | seq(b"le")
    | seq(b"in")
    | seq("≥".as_bytes())
    | seq("≤".as_bytes())
).convert(|arg0: &[u8]| String::from_utf8(arg0.to_vec()))
}

fn lparen<'a>() -> Parser<'a, u8, ()> {
    seq(b"(").discard()
}

fn rparen<'a>() -> Parser<'a, u8, ()> {
    seq(b")").discard()
}

fn real_number<'a>() -> Parser<'a, u8, f64> {
    let integer = one_of(b"123456789") - one_of(b"0123456789").repeat(0..) | sym(b'0');
    let frac = sym(b'.') + one_of(b"0123456789").repeat(1..);
    let exp = one_of(b"eE") + one_of(b"+-").opt() + one_of(b"0123456789").repeat(1..);
    let number = sym(b'-').opt() + integer + frac.opt() + exp.opt();
    number
        .collect()
        .convert(str::from_utf8)
        .convert(|s| f64::from_str(&s))
}

fn integer<'a>() -> Parser<'a, u8, u8> {
    one_of(b"123456789") - one_of(b"0123456789").repeat(0..) | sym(b'0')
}

fn str<'a>() -> Parser<'a, u8, String> {
    (sym(b'"') * none_of(b"\"").repeat(0..) - sym(b'"')).convert(String::from_utf8)
}

fn bool<'a>() -> Parser<'a, u8, Value> {
    ((seq(b"t") | seq(b"T"))
        + (seq(b"r") | seq(b"R"))
        + (seq(b"u") | seq(b"U"))
        + (seq(b"e") | seq(b"E")))
    .map(|_| Value::Bool(true))
        | ((seq(b"f") | seq(b"F"))
            + (seq(b"a") | seq(b"A"))
            + (seq(b"l") | seq(b"L"))
            + (seq(b"s") | seq(b"S"))
            + (seq(b"e") | seq(b"E")))
        .map(|_| Value::Bool(true))
}

fn none<'a>() -> Parser<'a, u8, u8> {
    ((seq(b"n") | seq(b"N"))
        + (seq(b"o") | seq(b"O"))
        + (seq(b"n") | seq(b"N"))
        + (seq(b"e") | seq(b"E")))
    .map(|_| 0)
}

fn value<'a>() -> Parser<'a, u8, Value> {
    space()
        * (real_number().map(|f| Value::from(f))
            | integer().map(|i| Value::Number(i.into()))
            | str().map(|s| Value::String(s))
            | bool()
            | none().map(|_| Value::Null)
            | ident().map(|s| Value::String(s))
        )
        - space()
}

fn and<'a>() -> Parser<'a, u8, u8> {
    ((seq(b"a") | seq(b"A")) + (seq(b"n") | seq(b"N")) + (seq(b"d") | seq(b"D"))).map(|_| 0)
}

fn or<'a>() -> Parser<'a, u8, u8> {
    ((seq(b"o") | seq(b"O")) + (seq(b"r") | seq(b"R"))).map(|_| 0)
}

fn and_or<'a>() -> Parser<'a, u8, String> {
    and().map(|_| "AND".into()) | or().map(|_| "OR".into())
}

fn property_val<'a>() -> Parser<'a, u8, Value> {
    space()
        * ((lparen() * list(value(), sym(b',').opt() * space()) - rparen())
            .map(|g| json!(g))
            | value())
        - space()
}

fn boolean_condition<'a>() -> Parser<'a, u8, Value> {
    space()
        * ((property_val() + operator() + property_val())
            .map(|((lval, op), rval)| json!({"name": lval, "op": op, "value": rval}))
            | (lparen() * call(boolean_expression) - rparen()).map(|boolean_expression| {
                json!(boolean_expression)
            }))
        - space()
}

fn and_expression<'a>() -> Parser<'a, u8, Value> {
    (boolean_condition() + and() * call(boolean_condition)).map(
        |(boolean_condition, and_expression)| json!({
            "AND": [boolean_condition, and_expression]
        })
    ) | boolean_condition()
}

fn or_expression<'a>() -> Parser<'a, u8, Value> {
    (and_expression() + or() * call(and_expression)).map(
        |(and_expression, or_expression)| json!({
            "OR": [and_expression, or_expression]
        })
    ) | and_expression()
}

fn boolean_expression<'a>() -> Parser<'a, u8, Value> {
    (boolean_condition() + and_or() + call(boolean_expression)).map(
        |((boolean_condition, and_or_initial), boolean_expression)| json!({
            and_or_initial: [boolean_condition, boolean_expression]
        })
    ) | boolean_condition()
}

// Define the main parser function
pub fn parse_logic_expr<'a>(input: &str) -> Result<Value, pom::Error> {
    //(space() * boolean_expression() - end()).parse(input.as_bytes())
    (space() * or_expression() - end()).parse(input.as_bytes())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_logic() -> Result<(), pom::Error> {

        let exprs = [
            (r#"foo = bar"#, r#"{"name":"foo","op":"=","value":"bar"}"#),
            (r#"foo = bar AND baz > 10"#, r#"{"AND":[{"name":"foo","op":"=","value":"bar"},{"name":"baz","op":">","value":10.0}]}"#),
            (r#"(foo == bar or baz > 10) AND baz <= 5"#, r#"{"AND":[{"OR":[{"name":"foo","op":"==","value":"bar"},{"name":"baz","op":">","value":10.0}]},{"name":"baz","op":"<=","value":5.0}]}"#),
            (r#"baz <= 5 AND (foo == bar or baz > 10)"#, r#"{"AND":[{"name":"baz","op":"<=","value":5.0},{"OR":[{"name":"foo","op":"==","value":"bar"},{"name":"baz","op":">","value":10.0}]}]}"#),
            (r#"foo == bar AND baz > 10 Or baz <= 5"#, r#"{"OR":[{"AND":[{"name":"foo","op":"==","value":"bar"},{"name":"baz","op":">","value":10.0}]},{"name":"baz","op":"<=","value":5.0}]}"#),
            (r#"baz <= 5 or foo == bar AND baz > 10  "#, r#"{"OR":[{"name":"baz","op":"<=","value":5.0},{"AND":[{"name":"foo","op":"==","value":"bar"},{"name":"baz","op":">","value":10.0}]}]}"#),
        ];
        for (expr, expected) in exprs.iter() {
            let result = parse_logic_expr(expr)?;
            assert_eq!(result, serde_json::from_str::<serde_json::Value>(expected).unwrap());
        }
        Ok(())
    }
}