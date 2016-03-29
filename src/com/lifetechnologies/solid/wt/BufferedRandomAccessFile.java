package com.lifetechnologies.solid.wt;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.io.File;

/**
 * User: tuchbb
 * Date: Oct 20, 2008
 * Time: 9:38:18 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 *
 * Credit is due to JavaWorld.com's "Java Tip 26" and Nick Zhang for the general design of this class. 
 */
public class BufferedRandomAccessFile extends RandomAccessFile {

    private int BUFFER_SIZE;

    byte buffer[];
    int indexBufferEnd = 0;
    int indexCurrentPositionInBuffer = 0;
    long positionInFile = 0;

    public BufferedRandomAccessFile(File file, String mode) throws IOException {
        super(file, mode);
        invalidate();
        this.BUFFER_SIZE = 1024 * 64;
        this.buffer = new byte[this.BUFFER_SIZE];
    }

    public BufferedRandomAccessFile(File file, String mode, int bufsize) throws IOException {
        super(file, mode);
        invalidate();
        this.BUFFER_SIZE = bufsize;
        this.buffer = new byte[this.BUFFER_SIZE];
    }


    public final int read() throws IOException {
        if (this.indexCurrentPositionInBuffer >= this.indexBufferEnd) {
            if (this.readBuffer() < 0)
                return -1;
        }
        if (this.indexBufferEnd == 0) {
            return -1;
        } else {
            return this.buffer[indexCurrentPositionInBuffer++];
        }
    }

    private int readBuffer() throws IOException {
        int numberOfBytesReadIntoBuffer = super.read(this.buffer, 0, this.BUFFER_SIZE);
        if (numberOfBytesReadIntoBuffer >= 0) {
            this.positionInFile += numberOfBytesReadIntoBuffer;
            this.indexBufferEnd = numberOfBytesReadIntoBuffer;
            this.indexCurrentPositionInBuffer = 0;
        }
        return numberOfBytesReadIntoBuffer;
    }

    private void invalidate() throws IOException {
        this.indexBufferEnd = 0;
        this.indexCurrentPositionInBuffer = 0;
        this.positionInFile = super.getFilePointer();
    }

    public int read(byte b[], int off, int len) throws IOException {
        int leftover = this.indexBufferEnd - this.indexCurrentPositionInBuffer;
        if (len <= leftover) {
            System.arraycopy(this.buffer, this.indexCurrentPositionInBuffer, b, off, len);
            this.indexCurrentPositionInBuffer += len;
            return len;
        }
        for (int i = 0; i < len; i++) {
            int c = this.read();
            if (c != -1)
                b[off + i] = (byte) c;
            else {
                return i;
            }
        }
        return len;
    }

    public long getFilePointer() throws IOException {
        long l = this.positionInFile;
        return (l - this.indexBufferEnd + this.indexCurrentPositionInBuffer);
    }

    public void seek(long pos) throws IOException {
        int n = (int) (positionInFile - pos);
        if (n >= 0 && n <= indexBufferEnd) {
            indexCurrentPositionInBuffer = indexBufferEnd - n;
        } else {
            super.seek(pos);
            invalidate();
        }
    }

    public final String readNextLine() throws IOException {
        String line = null;

        if (this.indexBufferEnd - this.indexCurrentPositionInBuffer <= 0) {
            if (readBuffer() < 0)
                return null;
                //throw new IOException("Error filling buffer.");

        }

        int indexOfLineEnd = -1;
        for (int i = this.indexCurrentPositionInBuffer; i < this.indexBufferEnd; i++) {
            if (buffer[i] == '\n') {
                indexOfLineEnd = i;
                break;
            }
        }

        if (indexOfLineEnd < 0) {
            StringBuffer input = new StringBuffer(256);
            int c;
            while (((c = read()) != -1) && (c != '\n')) {
                input.append((char) c);
            }
            if ((c == -1) && (input.length() == 0)) {
                return null;
            }
            return input.toString();
        }

        if (indexOfLineEnd > 0 && this.buffer[indexOfLineEnd - 1] == '\r')
            line = new String(this.buffer, this.indexCurrentPositionInBuffer, indexOfLineEnd - this.indexCurrentPositionInBuffer - 1);  // faster
            //line = new String(this.buffer, 0, this.indexCurrentPositionInBuffer, indexOfLineEnd - this.indexCurrentPositionInBuffer - 1);
        else
            line = new String(this.buffer, this.indexCurrentPositionInBuffer, indexOfLineEnd - this.indexCurrentPositionInBuffer);
            //line = new String(this.buffer, 0, this.indexCurrentPositionInBuffer, indexOfLineEnd - this.indexCurrentPositionInBuffer);

        this.indexCurrentPositionInBuffer = indexOfLineEnd + 1;

        return line;
    }


}
