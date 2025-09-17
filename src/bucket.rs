use std::io::{Read, Write};
use stream_vbyte::decode::decode;
use stream_vbyte::encode::encode;
use stream_vbyte::scalar::Scalar;

type Gid = u32;

#[derive(Debug, Clone, Default)]
pub struct Bucket {
    gids: Vec<Gid>,
}

impl Bucket {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn add_gid(&mut self, id: Gid) {
        self.gids.push(id);
    }

    pub fn get_vect_gid(&self) -> &Vec<Gid> {
        &self.gids
    }

    pub fn get_nb_gid(&self) -> u32 {
        self.gids.len() as u32
    }

    pub fn final_compression(&mut self) {
        // This is a no-op as the data is kept in a Vec.
        // This is kept for API compatibility with the builder pattern.
    }

    pub fn serialize<W: Write>(&mut self, writer: &mut W) -> std::io::Result<()> {
        let mut deltas = self.gids.clone();
        if !deltas.is_empty() {
            deltas.sort_unstable(); // Sort for delta encoding
            let mut prev = deltas[0];
            for i in 1..deltas.len() {
                let current = deltas[i];
                deltas[i] = current.wrapping_sub(prev);
                prev = current;
            }
        }
        
        // Allocate a buffer large enough for the worst-case scenario.
        let mut compressed_data = vec![0u8; deltas.len() * 5];
        let bytes_written = encode::<Scalar>(&deltas, &mut compressed_data);
        compressed_data.truncate(bytes_written);

        let num_elements = self.gids.len() as u32;
        let compressed_size = compressed_data.len() as u32;

        writer.write_all(&compressed_size.to_le_bytes())?;
        writer.write_all(&num_elements.to_le_bytes())?;
        writer.write_all(&compressed_data)?;
        Ok(())
    }

    pub fn deserialize<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut u32_buf = [0u8; 4];

        reader.read_exact(&mut u32_buf)?;
        let compressed_size = u32::from_le_bytes(u32_buf);

        reader.read_exact(&mut u32_buf)?;
        let num_elements = u32::from_le_bytes(u32_buf) as usize;

        let mut compressed_data = vec![0u8; compressed_size as usize];
        if compressed_size > 0 {
            reader.read_exact(&mut compressed_data)?;
        }
        
        let mut deltas = vec![0u32; num_elements];
        if num_elements > 0 {
            decode::<Scalar>(&compressed_data, num_elements, &mut deltas);
        }

        if !deltas.is_empty() {
            for i in 1..deltas.len() {
                deltas[i] = deltas[i].wrapping_add(deltas[i-1]);
            }
        }

        Ok(Self { gids: deltas })
    }
}
