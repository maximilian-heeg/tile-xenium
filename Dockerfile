FROM rust AS build
WORKDIR /app
COPY . .
RUN cargo build --release 
RUN strip target/release/tile-xenium
ENTRYPOINT ["/app/target/release/tile-xenium"]
CMD ["--help"]


FROM debian:bookworm-slim
WORKDIR /app
COPY --from=build /app/target/release/tile-xenium /app/
ENV PATH="${PATH}:/app/"
ENTRYPOINT ["/app/tile-xenium"]
CMD ["--help"]