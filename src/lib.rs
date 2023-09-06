use clap::Parser;
use polars::prelude::*;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::prelude::*;

/// Create tiles, filter transcripts from transcripts.csv based
/// on Q-Score threshold. Remove negative controls
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// The path to the transcripts.csv or .parquet file produced
    /// by Xenium.
    in_file: String,

    /// The minimum Q-Score to pass filtering.
    #[arg(long, default_value_t = 20.0)]
    min_qv: f64,

    /// The width of the tiles
    #[arg(long, default_value_t = 4000.0)]
    width: f64,

    /// The height of the tiles
    #[arg(long, default_value_t = 4000.0)]
    height: f64,

    /// Overlap betweent the tiles
    #[arg(long, default_value_t = 500.0)]
    overlap: f64,

    /// Minimal number of transcripts per tile. If number is less, the tile will be expanded by overlap in all directions.
    #[arg(long, default_value_t = 100000)]
    minimal_transcripts: usize,

    /// Only keep transcripts that are in the nucelus.
    /// All other transcripts will not be assigned to a cell.
    #[arg(long, default_value_t = false)]
    nucleus_only: bool,

    /// Path for outout file. The output file will be names after the following schema:  
    ///   X{x-min}-{x-max}_Y{y-min}-{y-max}_filtered_transcripts_nucleus_only_{nucleus_only}.csv  
    ///   E.g.: X0-24000_Y0-24000_filtered_transcripts_nucleus_only_false.csv
    #[arg(long, default_value = ".", verbatim_doc_comment)]
    out_dir: String,
}

#[derive(Debug)]
struct MyError(String);

impl fmt::Display for MyError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Error: {}", self.0)
    }
}
impl Error for MyError {}

pub fn run(args: Args) -> Result<(), Box<dyn Error>> {
    if args.width <= args.overlap {
        return Err(Box::new(MyError(
            "The width of a tile cannot be smaller than the overlap.".into(),
        )));
    }
    if args.width <= args.overlap {
        return Err(Box::new(MyError(
            "The height of a tile cannot be smaller than the overlap.".into(),
        )));
    }

    println!("Reading file");
    let mut transcripts: LazyFrame = if args.in_file.ends_with(".parquet") {
        LazyFrame::scan_parquet(&args.in_file, Default::default())?
    } else if args.in_file.ends_with(".csv") {
        LazyCsvReader::new(&args.in_file)
            .has_header(true)
            .finish()?
    } else {
        return Err(Box::new(MyError(
            "Input file should be either CSV or Parquet.".into(),
        )));
    };

    transcripts = transcripts.select([
        col("transcript_id"),
        col("cell_id"),
        col("overlaps_nucleus"),
        col("feature_name"),
        col("x_location"),
        col("y_location"),
        col("z_location"),
        col("qv"),
    ]);

    // Decode the binary data column
    transcripts = transcripts.with_columns([
        col("feature_name").cast(DataType::Utf8),
        col("cell_id").cast(DataType::Utf8),
    ]);

    // Remove the non-gene features
    println!("Remove non-gene features");
    transcripts = filter_transcripts(transcripts);
    let total_transcripts = get_row_count(transcripts.clone());
    if total_transcripts < args.minimal_transcripts.try_into()? {
        return Err(Box::new(MyError(
            "The number of transcripts that remain after excluding \
            the non-gene transcripts is lower than the required minimal \
            transciprt number per tile. Please consider adjusting that value."
                .into(),
        )));
    }

    println!("Decoding cell ids");
    // Change cell ids from UNASSINGED to 0
    transcripts = transcripts.with_columns([when(col("cell_id").eq(lit("UNASSIGNED")))
        .then(lit("0"))
        .otherwise(col("cell_id"))
        .alias("cell_id")]);

    // remove cell assigment for transcripts, that are not in the nucleus
    if args.nucleus_only {
        transcripts = transcripts.with_columns([when(col("overlaps_nucleus").eq(lit(0)))
            .then(lit("0"))
            .otherwise(col("cell_id"))
            .alias("cell_id")]);
    }
    // Decode the 10x cell id into an integer
    transcripts = decode_cells(transcripts);
    let df = transcripts.collect()?;
    println!("{}", df);

    println!("Get limits");
    let (x_min, x_max, y_min, y_max) = get_limits(&df);
    println!("... x: {x_min} - {x_max}");
    println!("... y: {y_min} - {y_max}");
    // Create tiltes
    println!("Create tiles");
    std::fs::create_dir_all(&args.out_dir)?;
    let mut y = y_min;
    let mut x = x_min;
    while y <= y_max {
        while x <= x_max {
            let start_x = x;
            let end_x = x + args.width;
            let start_y = y;
            let end_y = y + args.height;

            println!(
                "... Trying to create tile from X= {start_x} - {end_x} and Y= {start_y} - {end_y}"
            );
            let file = write_tile(
                df.clone().lazy(),
                start_x,
                end_x,
                start_y,
                end_y,
                &args.out_dir,
                args.minimal_transcripts,
                args.overlap,
            );
            println!("... ... Tile created: {file}");

            x = x + args.width - args.overlap;
        }
        x = x_min;
        y = y + args.height - args.overlap;
    }

    // Dump parameters
    let outfile = format!("{}/params.txt", args.out_dir);
    let mut file = File::create(outfile)?;
    write!(file, "{:?}", args)?;

    Ok(())
}

/// Write tiles
fn write_tile(
    df: LazyFrame,
    start_x: f64,
    end_x: f64,
    start_y: f64,
    end_y: f64,
    outdir: &String,
    minimal_transcripts: usize,
    increment: f64,
) -> String {
    let filtered_df = df
        .clone()
        .filter(col("x_location").gt_eq(start_x))
        .filter(col("x_location").lt_eq(end_x))
        .filter(col("y_location").gt_eq(start_y))
        .filter(col("y_location").lt_eq(end_y));

    let count = get_row_count(filtered_df.clone());

    if count < minimal_transcripts as u32 {
        print!("... ... Not enought transcripts. Adding buffer to tile.\n");
        return write_tile(
            df,
            start_x - increment,
            end_x + increment,
            start_y - increment,
            end_y + increment,
            outdir,
            minimal_transcripts,
            increment,
        );
    }

    // Create output file
    let mut df = filtered_df.collect().unwrap();
    let outfile = format!(
        "{}/X{}-{}_Y{}-{}_filtered_transcripts.csv",
        outdir, start_x, end_x, start_y, end_y
    );
    let mut output_file: File = File::create(&outfile).unwrap();
    CsvWriter::new(&mut output_file).finish(&mut df).unwrap();
    return outfile;
}

/// Get the minimal and maximal values of x_location and y_location
/// values are rounded to the smallest integer number greater than the value for the maximal values
/// values are rounded to the biggest integer number lower than the value for the minimal values
fn get_limits(df: &DataFrame) -> (f64, f64, f64, f64) {
    let x_min: f64 = df.column("x_location").unwrap().min().unwrap();

    let x_max: f64 = df.column("x_location").unwrap().max().unwrap();

    let y_min: f64 = df.column("y_location").unwrap().min().unwrap();

    let y_max: f64 = df.column("y_location").unwrap().max().unwrap();

    (x_min.floor(), x_max.ceil(), y_min.floor(), y_max.ceil())
}

fn get_row_count(df: LazyFrame) -> u32 {
    df.select([count().alias("count")])
        .collect()
        .unwrap()
        .column("count")
        .unwrap()
        .u32()
        .unwrap()
        .get(0)
        .unwrap()
}
/// Exclude control probes from the transcripts table.
fn filter_transcripts(t: LazyFrame) -> LazyFrame {
    t.filter(not(col("feature_name")
        .str()
        .starts_with(lit("NegControlProbe_"))))
        .filter(not(col("feature_name")
            .str()
            .starts_with(lit("antisense_"))))
        .filter(not(col("feature_name")
            .str()
            .starts_with(lit("NegControlCodeword_"))))
        .filter(not(col("feature_name").str().starts_with(lit("BLANK_"))))
}

/// Decode the cell_id column
/// Set unassinged to zero
/// Set cell to zero if not in nucleus and option enabled
fn decode_cells(df: LazyFrame) -> LazyFrame {
    let result = df.with_columns([col("cell_id").map(
        |x| {
            let out: Series = x
                .iter()
                .map(|v: AnyValue<'_>| match v {
                    AnyValue::Utf8(val) => decode_cell_id(val),
                    _v => 0,
                })
                .collect();
            Ok(Some(out))
        },
        GetOutput::from_type(DataType::Utf8),
    )]);
    result
}

fn decode_cell_id(cell_id: &str) -> u32 {
    // Check if id can be converted to integer
    if let Ok(numeric_id) = cell_id.parse::<u32>() {
        return numeric_id;
    }
    let mut parts = cell_id.split("-");
    let shifted_hex_digits = parts.next().unwrap();
    let _dataset_suffix: usize = parts.next().unwrap().parse().unwrap();

    let hex_array = shifted_hex_to_hex_array(shifted_hex_digits);

    let integer_value = convert_hex_array_to_int(&hex_array);

    integer_value
}

fn shifted_hex_to_hex_array(shifted_hex_digits: &str) -> Vec<i32> {
    let hex_digits: Vec<i32> = shifted_hex_digits
        .chars()
        .map(|c| (c as i32) - ('a' as i32))
        .collect();

    hex_digits
}

fn convert_hex_array_to_int(hex_array: &[i32]) -> u32 {
    let hex_string: String = hex_array.iter().map(|x| format!("{:X}", x)).collect();

    let integer_value = u32::from_str_radix(&hex_string, 16).unwrap();
    integer_value
}

#[cfg(test)]
mod tests {
    use crate::filter_transcripts;

    use super::*;

    #[test]
    fn reverse_cell_id() {
        let cell_id = "ffkpbaba-1";

        let decoded = decode_cell_id(cell_id);
        assert_eq!(decoded, 1437536272);
    }

    #[test]
    fn check_decode_cells() {
        let df = df! [
            "cell_id" => ["0", "ffkpbaba-1", "ffkpbaba-1"],
            "overlaps_nucleus"     => [1,1,0]
        ]
        .unwrap()
        .lazy();

        let solution = df! [
            "cell_id" => [0, 1437536272, 1437536272],
            "overlaps_nucleus"     => [1,1,0]
        ]
        .unwrap();

        let res = decode_cells(df).collect().unwrap();
        assert_eq!(res, solution);
    }

    #[test]
    fn check_filter_transcripts() {
        let df = df! [
            "feature_name" => ["NegControlProbe_", "antisense_", "NegControlCodeword_", "BLANK_", "Gene"],
            "Other_col"         => [false, false, true, false, true]
        ].unwrap().lazy();

        let result = filter_transcripts(df).collect().unwrap();

        let solution = df! [
            "feature_name" => ["Gene"],
            "Other_col"         => [true]
        ]
        .unwrap();

        assert_eq!(result, solution);
    }

    #[test]
    fn check_limits() {
        let df = df! [
            "x_location" => [1.0, 2.0, 3.0, 2.5],
            "y_location" => [9.0, 8.7, 8.8, 8.9]
        ]
        .unwrap();

        assert_eq!(get_limits(&df), (1.0, 3.0, 8.0, 9.0))
    }
}
