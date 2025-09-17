use std::io::{Read, Write};
use stream_vbyte::decode::decode;
use stream_vbyte::encode::encode;
use stream_vbyte::scalar::Scalar;

#[derive(Debug, Clone, Default)]
pub struct BucketPos {
    positions: Vec<u16>,
}

impl BucketPos {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn add_pos(&mut self, pos: u16) {
        self.positions.push(pos);
    }

    pub fn get_vect_pos(&self) -> &Vec<u16> {
        &self.positions
    }

    pub fn get_nb_pos(&self) -> u32 {
        self.positions.len() as u32
    }

    pub fn final_compression(&mut self) {
        // No-op, kept for API compatibility.
    }

    pub fn serialize<W: Write>(&mut self, writer: &mut W) -> std::io::Result<()> {
        let data_u32: Vec<u32> = self.positions.iter().map(|&p| p as u32).collect();

        // Allocate a buffer large enough for the worst-case scenario.
        let mut compressed_data = vec![0u8; data_u32.len() * 5];
        let bytes_written = encode::<Scalar>(&data_u32, &mut compressed_data);
        compressed_data.truncate(bytes_written);

        let num_elements = self.positions.len() as u32;
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
        
        let mut positions_u32 = vec![0u32; num_elements];
        if num_elements > 0 {
            // Specify the Scalar decoder implementation
            decode::<Scalar>(&compressed_data, num_elements, &mut positions_u32);
        }

        let positions = positions_u32.into_iter().map(|p| p as u16).collect();
        Ok(Self { positions })
    }
}
