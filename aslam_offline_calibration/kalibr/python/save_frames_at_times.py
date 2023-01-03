"""
save frames from a rosbag at times to another rosbag
"""
import rosbag
import rospy
import sys


def main():
    if len(sys.argv) < 5:
        print("Usage: {} <bagfile> <imagetopic> <timefile> <outputbag>".format(sys.argv[0]))
        return

    inbag = sys.argv[1]
    topic = sys.argv[2]
    timefile = sys.argv[3]
    outbag = sys.argv[4]

    id2times = []
    with open(timefile) as s:
        for l in s:
            fields = l.split()
            id2times.append((int(fields[0]), rospy.Time(float(fields[1]))))

    # select and save the static images
    reader = rosbag.Bag(inbag, 'r')
    writer = rosbag.Bag(outbag, 'w')
    count = 0
    numframes = 0
    j = 0
    for _, msg, t in reader.read_messages(topics=[topic]):
        if (msg.header.stamp - id2times[j][1]).nsecs < 1000:
            writer.write(topic, msg, msg.header.stamp)
            count += 1
            j += 1
        if j >= len(id2times):
            break
        numframes += 1
    print('Saved {} frames out of {} frames to {}.'.format(count, numframes, outbag))
    if count != len(id2times):
        print("Warning: #required frames {} is more than the found {}.".format(len(id2times), count))
    writer.close()


if __name__ == '__main__':
    main()
